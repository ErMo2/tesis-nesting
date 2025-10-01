package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class GraspPlacer {

    // ===================== Parámetros GRASP =====================
    static final int    ITER_MAX            = 120;   // número máximo de iteraciones GRASP
    static final int    RCL_SIZE            = 10;    // tamaño de la lista restringida de candidatos
    static final long   MAX_MILLIS_TOTAL    = 45_000; // tiempo total permitido por corrida
    static final long   MAX_MILLIS_LOCAL    = 12_000; // tiempo para búsqueda local
    static final long   DEFAULT_SEED        = 42L;   // semilla default

    // parámetros de anclas / chequeos
    static final int    GRID_ANCHORS_X      = 10;
    static final int    GRID_ANCHORS_Y      = 6;
    static final int    MAX_ANCHORS_PER_TRY = 220;
    static final int    MAX_CHECKS_PER_PART = 1200;
    static final double NUDGE_STEP          = 0.35;
    static final int    NUDGE_STEPS_MAX     = 6;

    // scoring local
    static final double W_BEST_Y            = 1.0;
    static final double W_BEST_X            = 0.15;
    static final double W_BBOX_INFLATE      = 0.001;

    // búsqueda local
    static final int    REINSERT_TRIES      = 2;
    static final int    TWEAK_ANGLE_STEP    = 3;

    // lattice mode
    static final double GRID_CV_MAX         = 0.02;
    static final double EPS_EDGE            = 1e-9;

    // logging
    static final String LOG_FILE = "LogGrasp.log";
    static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");
    static final DecimalFormat DF = new DecimalFormat("#,##0.###");
    private static PrintStream tee;

    // ======= Tipos =======

    // representa una copia de una pieza
    static class PartInstance {
        final PartSpec spec; final int copyIdx;
        PartInstance(PartSpec s, int k) { spec = s; copyIdx = k; }
        String key() { return spec.id + "#" + copyIdx; }
    }

    // representa un punto de anclaje
    static class Anchor { final double x, y; Anchor(double x, double y){ this.x=x; this.y=y; } }

    // candidato durante la construcción (geom orientada + datos de colocación)
    static class Candidate {
        Geometry clearGeom;   // con clearance
        Geometry placedGeom;  // pieza real
        double   dx, dy;
        int      angleDeg;
        double   scoreVal;
    }

    // resultado de una corrida GRASP
    static class RunResult {
        List<PlacedPart> placed = Collections.emptyList();
        double bestArea;
        double bestY;
        int totalPieces;
        long millis;
        double areaSheet;
        double wasted;
        double utilizationPct;
        long seedUsed;
    }

    // ======= API =======

    // API principal: corre GRASP con semilla por defecto
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        RunResult rr = runOnce(sheetInset, specs, DEFAULT_SEED, true);
        return rr.placed;
    }

    // API con semilla personalizada
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs, long seed) {
        RunResult rr = runOnce(sheetInset, specs, seed, true);
        return rr.placed;
    }

    // API multi-start con seeds distintos
    public static List<PlacedPart> placeMulti(Polygon sheetInset, List<PartSpec> specs, int numRuns) {
        return placeMulti(sheetInset, specs, numRuns, System.nanoTime());
    }

    // API multi-start reproducible (baseSeed+i)
    public static List<PlacedPart> placeMulti(Polygon sheetInset, List<PartSpec> specs, int numRuns, long baseSeed) {
        if (numRuns <= 0) numRuns = 1;

        initLogger();
        long t0 = System.currentTimeMillis();
        log("==== GRASP Multi-Start RUNS=" + numRuns + " ====");

        RunResult best = null;
        for (int i=0; i<numRuns; i++) {
            long seed = baseSeed + i;
            RunResult rr = runOnce(sheetInset, specs, seed, false);

            log(String.format("RUN %d | seed=%d | colocadas=%d/%d | area=%s | util=%s%% | yMax=%s | dur=%d ms",
                    i+1, seed, rr.placed.size(), rr.totalPieces,
                    DF.format(rr.bestArea), DF.format(rr.utilizationPct), DF.format(rr.bestY), rr.millis));

            if (best == null ||
                rr.bestArea > best.bestArea ||
               (Math.abs(rr.bestArea - best.bestArea) < 1e-6 && rr.bestY < best.bestY)) {
                best = rr;
            }
        }

        long totalMs = System.currentTimeMillis() - t0;
        if (best == null) {
            log("No se pudo construir ninguna solución. Fin.");
            closeLogger();
            return Collections.emptyList();
        }

        log("------------------------------");
        log("RESUMEN GLOBAL (Multi-Start)");
        log("Mejor seed           : " + best.seedUsed);
        log("Piezas a posicionar  : " + best.totalPieces);
        log("Piezas colocadas     : " + best.placed.size());
        log("% aprovechamiento    : " + DF.format(best.utilizationPct) + " %");
        log("Área colocada        : " + DF.format(best.bestArea));
        log("Área plancha         : " + DF.format(best.areaSheet));
        log("Área desperdiciada   : " + DF.format(best.wasted));
        log("Ymax                 : " + DF.format(best.bestY));
        log("Tiempo total (multi) : " + totalMs + " ms (" + DF.format(totalMs/1000.0) + " s)");
        log("------------------------------");

        closeLogger();
        return best.placed;
    }
    // ======= Núcleo de una corrida GRASP =======
    // se encarga de expandir copias, construir soluciones y aplicar búsqueda local
    private static RunResult runOnce(Polygon sheetInset, List<PartSpec> specs, long seed, boolean handleLogger) {
        if (handleLogger) initLogger();
        final long t0 = System.currentTimeMillis();
        if (handleLogger) {
            log("==== GraspPlacer RUN (single) " + new Date() + " ====");
            log("INI | PartSpec: " + specs.size() + " | seed=" + seed);
        }

        // expandir las copias de cada pieza
        List<PartInstance> instances = expandInstances(specs);
        if (instances.isEmpty()) {
            if (handleLogger) { log("No hay copias a colocar. Fin."); closeLogger(); }
            RunResult empty = new RunResult();
            empty.placed = Collections.emptyList();
            empty.bestArea = 0;
            empty.bestY = 0;
            empty.totalPieces = 0;
            empty.millis = System.currentTimeMillis() - t0;
            empty.areaSheet = sheetInset.getArea();
            empty.wasted = empty.areaSheet;
            empty.utilizationPct = 0;
            empty.seedUsed = seed;
            return empty;
        }

        // calcular opciones de ángulo por pieza
        Map<PartSpec,int[]> angleOptions = computeAngleOptions(specs);

        // loop GRASP: varias iteraciones buscando la mejor solución
        Random rng = new Random(seed);
        List<PlacedPart> bestSolution = Collections.emptyList();
        double bestArea = -1, bestY = Double.POSITIVE_INFINITY;

        for (int it = 0; it < ITER_MAX; it++) {
            if (System.currentTimeMillis() - t0 > MAX_MILLIS_TOTAL) { 
                if (handleLogger) log("Tiempo agotado en iter=" + it); 
                break; 
            }
            long tIt = System.currentTimeMillis();

            // ordenar piezas: grandes primero, pero con mezcla aleatoria por bloques
            List<PartInstance> order = new ArrayList<>(instances);
            order.sort((a,b) -> Double.compare(b.spec.poly.getArea(), a.spec.poly.getArea()));
            blockShuffle(order, rng, 6);

            // construcción codiciosa aleatorizada
            DecodeResult dr = greedyRandomDecode(sheetInset, order, angleOptions, rng);

            // búsqueda local si se colocó algo y queda tiempo
            long tLocalBudget = Math.max(0, MAX_MILLIS_LOCAL - (System.currentTimeMillis()-tIt));
            if (!dr.placed.isEmpty() && tLocalBudget > 0) {
                localSearch(sheetInset, dr, angleOptions, rng, tLocalBudget);
            }

            // métricas de la solución actual
            double itArea = totalArea(dr.placed);
            double itYmax = usedYmax(dr.placed);

            // si mejora al mejor hasta ahora, lo guardo
            if (itArea > bestArea || (Math.abs(itArea-bestArea) < 1e-6 && itYmax < bestY)) {
                bestArea = itArea; 
                bestY = itYmax; 
                bestSolution = dr.placed;
            }

            if (handleLogger) {
                log(String.format("IT %d | colocadas=%d | area=%s | yMax=%s | dur=%d ms",
                        it, dr.placed.size(), DF.format(itArea), DF.format(itYmax),
                        (System.currentTimeMillis()-tIt)));
            }
        }

        // ====== RESUMEN DE LA CORRIDA ======
        int totalPiezas        = instances.size();
        int colocadas          = bestSolution.size();
        double areaPlancha     = sheetInset.getArea();
        double areaColocada    = totalArea(bestSolution);
        double areaDesperdiciada = Math.max(0.0, areaPlancha - areaColocada);
        double aprovechamiento = areaPlancha > 0 ? (100.0 * areaColocada / areaPlancha) : 0.0;
        long durTotalMs        = System.currentTimeMillis() - t0;

        if (handleLogger) {
            log("------------------------------");
            log("RESUMEN FINAL (GRASP - single run)");
            log("Seed               : " + seed);
            log("Piezas a posicionar: " + totalPiezas);
            log("Piezas colocadas   : " + colocadas);
            log("% aprovechamiento  : " + DF.format(aprovechamiento) + " %");
            log("Área colocada      : " + DF.format(areaColocada));
            log("Área plancha       : " + DF.format(areaPlancha));
            log("Área desperdiciada : " + DF.format(areaDesperdiciada));
            log("Ymax               : " + DF.format(bestY));
            log("Tiempo total       : " + durTotalMs + " ms (" + DF.format(durTotalMs/1000.0) + " s)");
            log("------------------------------");
            log("BEST | colocadas=" + colocadas + " | area=" + DF.format(areaColocada) + " | yMax=" + DF.format(bestY));
            closeLogger();
        }

        // armar objeto con los resultados
        RunResult rr = new RunResult();
        rr.placed = bestSolution;
        rr.bestArea = areaColocada;
        rr.bestY = bestY;
        rr.totalPieces = totalPiezas;
        rr.millis = durTotalMs;
        rr.areaSheet = areaPlancha;
        rr.wasted = areaDesperdiciada;
        rr.utilizationPct = aprovechamiento;
        rr.seedUsed = seed;
        return rr;
    }
    // ======= Estructura para guardar la decodificación =======
    // contiene la lista de piezas colocadas, geometrías con clearance y la bbox global
    static class DecodeResult {
        final List<PlacedPart> placed = new ArrayList<>();
        final List<Geometry>   placedClear = new ArrayList<>();
        double globMinX, globMinY, globMaxX, globMaxY;
    }

    // ======= Construcción codiciosa aleatorizada (GRASP) =======
    // genera una solución parcial/total probando ángulos y anclas con RCL
    static DecodeResult greedyRandomDecode(Polygon sheetInset,
                                           List<PartInstance> order,
                                           Map<PartSpec,int[]> angleOptions,
                                           Random rng)
    {
        DecodeResult res = new DecodeResult();
        Envelope env = sheetInset.getEnvelopeInternal();
        double minX = env.getMinX(), minY = env.getMinY(), maxX = env.getMaxX(), maxY = env.getMaxY();
        res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;

        // -------- Detectar LATTICE_MODE --------
        boolean GRID_MODE = false;
        double cellW = 0, cellH = 0;
        {
            double rectArea = (maxX - minX) * (maxY - minY);
            boolean sheetRectangular = Math.abs(sheetInset.getArea() - rectArea) < 1e-6;

            int limit = Math.min(50, order.size());
            if (sheetRectangular && limit > 0) {
                double[] ws = new double[limit], hs = new double[limit];
                for (int i=0;i<limit;i++){
                    PartSpec sp = order.get(i).spec;
                    int[] allowed = angleOptions.getOrDefault(sp, new int[]{0});
                    int angDeg = allowed.length>0 ? allowed[0] : 0;
                    Geometry g0  = GeomUtils.orient(sp.poly, angDeg, false);
                    Geometry gc0 = (sp.clearance > 0) ? g0.buffer(sp.clearance, 8) : g0;
                    Envelope e   = gc0.getEnvelopeInternal();
                    ws[i] = e.getWidth(); hs[i] = e.getHeight();
                }
                if (cv(ws) <= GRID_CV_MAX && cv(hs) <= GRID_CV_MAX) {
                    GRID_MODE = true; 
                    cellW = median(ws); 
                    cellH = median(hs);
                }
            }
        }

        // si se detecta modo rejilla, usar fast-path
        if (GRID_MODE) {
            log("  GRID fast-path activo (GRASP) | cellW="+DF.format(cellW)+" | cellH="+DF.format(cellH));
            gridFillRaster(sheetInset, res.placed, res.placedClear, order, angleOptions,
                           minX, minY, maxX, maxY, cellW, cellH);

            // actualizar bbox global según lo colocado
            res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;
            for (Geometry cc : res.placedClear) {
                Envelope ce = cc.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());
            }
            return res; // <<< en GRID_MODE no seguimos con heurística normal
        }

        // -------- Flujo normal (NO rejilla) ----------
        ArrayList<Anchor> anchors = initialAnchors(minX, minY, maxX, maxY);

        int pos = 0;
        for (PartInstance inst : order) {
            PartSpec sp = inst.spec;
            int[] angles = angleOptions.getOrDefault(sp, new int[]{0});
            List<Candidate> rcl = new ArrayList<>(RCL_SIZE);
            int checks = 0;

            // muestreo de anclas para esta pieza
            List<Anchor> sampleAnchors = sampleAnchors(anchors, rng, MAX_ANCHORS_PER_TRY);

            // probar cada ángulo permitido
            for (int angDeg : angles) {
                Geometry gBase = GeomUtils.orient(sp.poly, angDeg, false);
                Geometry gClr  = sp.clearance > 0 ? gBase.buffer(sp.clearance, 8) : gBase;

                Envelope pe = gClr.getEnvelopeInternal();
                for (Anchor a : sampleAnchors) {
                    if (++checks > MAX_CHECKS_PER_PART) break;

                    double dxa = a.x - pe.getMinX();
                    double dya = a.y - pe.getMinY();

                    double txMin = pe.getMinX() + dxa, txMax = pe.getMaxX() + dxa;
                    double tyMin = pe.getMinY() + dya, tyMax = pe.getMaxY() + dya;
                    if (txMin < minX-EPS_EDGE || txMax > maxX+EPS_EDGE || 
                        tyMin < minY-EPS_EDGE || tyMax > maxY+EPS_EDGE) continue;

                    Geometry cc = GeomUtils.translate(gClr, dxa, dya);
                    if (!sheetInset.covers(cc)) continue; // se permite contacto borde-borde

                    boolean inter = false;
                    for (Geometry pg : res.placedClear) { 
                        if (collideStrict(pg, cc)) { inter = true; break; } 
                    }
                    if (inter) continue;

                    // mini compactación hacia la izquierda
                    double[] nudged = tryNudge(sheetInset, cc, res.placedClear, -NUDGE_STEP, 0, NUDGE_STEPS_MAX);
                    cc = GeomUtils.translate(cc, nudged[0], nudged[1]);
                    Geometry placed = GeomUtils.translate(gBase, dxa+nudged[0], dya+nudged[1]);

                    // calcular score local
                    double val = localScore(res, cc);

                    Candidate cand = new Candidate();
                    cand.clearGeom  = cc;
                    cand.placedGeom = placed;
                    cand.dx = dxa + nudged[0];
                    cand.dy = dya + nudged[1];
                    cand.angleDeg = angDeg;
                    cand.scoreVal = val;

                    insertIntoRCL(rcl, cand, RCL_SIZE);
                }
            }

            // si hay candidatos, elegir uno aleatorio de la RCL
            if (!rcl.isEmpty()) {
                Candidate pick = rcl.get(rng.nextInt(rcl.size()));

                PlacedPart pp = new PlacedPart();
                pp.id = sp.id;
                pp.polyPlaced = pick.placedGeom;
                pp.angleDeg = pick.angleDeg;
                pp.mirrored = false;
                pp.dx = pick.dx; pp.dy = pick.dy;

                res.placed.add(pp);
                res.placedClear.add(pick.clearGeom);

                // actualizar bbox global
                Envelope ce = pick.clearGeom.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());

                // generar nuevas anclas locales a partir de la pieza
                addLocalAnchors(anchors, ce, minX, minY, maxX, maxY);
            }

            if ((++pos % 10) == 0) {
                log("  · placed=" + res.placed.size() + " / " + order.size() + " | anchors=" + anchors.size());
            }
        }

        return res;
    }
    // ======= Búsqueda local =======
    // intenta mejorar la solución reinsertando piezas con tweaks de ángulo
    static void localSearch(Polygon sheetInset,
                            DecodeResult dr,
                            Map<PartSpec,int[]> angleOptions,
                            Random rng,
                            long millisBudget)
    {
        long t0 = System.currentTimeMillis();
        if (dr.placed.size() <= 1) return;

        // se hacen varios intentos de reinserción
        for (int pass=0; pass<REINSERT_TRIES; pass++) {
            // recorremos las piezas al revés para ir quitando y probando
            for (int i = dr.placed.size()-1; i >= 0; i--) {
                if (System.currentTimeMillis()-t0 > millisBudget) return;

                // quitar pieza candidata (victim)
                PlacedPart victim = dr.placed.remove(i);
                Geometry victimClear = dr.placedClear.remove(i);

                // construir anclas: plancha + piezas ya colocadas
                Envelope env = sheetInset.getEnvelopeInternal();
                ArrayList<Anchor> anchors = initialAnchors(env.getMinX(), env.getMinY(), env.getMaxX(), env.getMaxY());
                for (Geometry gc : dr.placedClear) {
                    addLocalAnchors(anchors, gc.getEnvelopeInternal(), env.getMinX(), env.getMinY(), env.getMaxX(), env.getMaxY());
                }

                // buscar la spec original de la pieza
                PartSpec sp = findSpec(angleOptions, victim.id);
                int[] baseAngles = angleOptions.getOrDefault(sp, new int[]{0});
                int[] angles = maybeAugmentAngles(baseAngles, sp);

                Candidate best = null;
                List<Anchor> sample = sampleAnchors(anchors, rng, MAX_ANCHORS_PER_TRY);

                // probar posibles ángulos de reinserción
                for (int angDeg : angles) {
                    Geometry gBase = GeomUtils.orient(sp.poly, angDeg, false);
                    Geometry gClr  = sp.clearance > 0 ? gBase.buffer(sp.clearance, 8) : gBase;
                    Envelope pe = gClr.getEnvelopeInternal();

                    int checks = 0;
                    for (Anchor a : sample) {
                        if (++checks > MAX_CHECKS_PER_PART) break;

                        double dxa = a.x - pe.getMinX();
                        double dya = a.y - pe.getMinY();

                        double txMin = pe.getMinX() + dxa, txMax = pe.getMaxX() + dxa;
                        double tyMin = pe.getMinY() + dya, tyMax = pe.getMaxY() + dya;
                        if (txMin < env.getMinX()-EPS_EDGE || txMax > env.getMaxX()+EPS_EDGE ||
                            tyMin < env.getMinY()-EPS_EDGE || tyMax > env.getMaxY()+EPS_EDGE) continue;

                        Geometry cc = GeomUtils.translate(gClr, dxa, dya);
                        if (!sheetInset.covers(cc)) continue; // contacto permitido en borde

                        boolean inter = false;
                        for (Geometry pg : dr.placedClear) { 
                            if (collideStrict(pg, cc)) { inter = true; break; } 
                        }
                        if (inter) continue;

                        // mini compactación hacia abajo/izquierda
                        double[] nudged = tryNudge(sheetInset, cc, dr.placedClear, -NUDGE_STEP, 0, NUDGE_STEPS_MAX);
                        cc = GeomUtils.translate(cc, nudged[0], nudged[1]);
                        Geometry placed = GeomUtils.translate(gBase, dxa+nudged[0], dya+nudged[1]);

                        // score local enfocado en la bbox
                        double val = localScoreBBoxOnly(sheetInset, dr, cc);

                        // actualizar mejor candidato
                        if (best == null || val < best.scoreVal) {
                            best = new Candidate();
                            best.clearGeom = cc; 
                            best.placedGeom = placed;
                            best.dx = dxa+nudged[0]; 
                            best.dy = dya+nudged[1];
                            best.angleDeg = angDeg; 
                            best.scoreVal = val;
                        }
                    }
                }

                // decidir si se acepta la reinserción
                if (best != null) {
                    double oldY = usedYmax(dr.placed); 
                    double oldArea = totalArea(dr.placed);
                    double newY, newArea;
                    {
                        // simular insertar el nuevo
                        PlacedPart pp = new PlacedPart();
                        pp.id = sp.id; 
                        pp.polyPlaced = best.placedGeom; 
                        pp.angleDeg = best.angleDeg; 
                        pp.mirrored = false; 
                        pp.dx = best.dx; 
                        pp.dy = best.dy;

                        dr.placed.add(pp); 
                        dr.placedClear.add(best.clearGeom);
                        newY = usedYmax(dr.placed); 
                        newArea = totalArea(dr.placed);

                        // revertir simulación
                        dr.placed.remove(dr.placed.size()-1); 
                        dr.placedClear.remove(dr.placedClear.size()-1);
                    }

                    // aceptar si no empeora
                    if (newY <= oldY + 1e-6 || newArea >= oldArea - 1e-6) {
                        PlacedPart pp = new PlacedPart();
                        pp.id = sp.id; 
                        pp.polyPlaced = best.placedGeom; 
                        pp.angleDeg = best.angleDeg; 
                        pp.mirrored = false; 
                        pp.dx = best.dx; 
                        pp.dy = best.dy;
                        dr.placed.add(i, pp); 
                        dr.placedClear.add(i, best.clearGeom);
                    } else {
                        // si no mejora, se devuelve la pieza original
                        dr.placed.add(i, victim); 
                        dr.placedClear.add(i, victimClear);
                    }
                } else {
                    // si no hubo mejor candidato, devolver la pieza original
                    dr.placed.add(i, victim); 
                    dr.placedClear.add(i, victimClear);
                }
            }
        }
    }
    // ======= Helpers / utilidades =======

    // inicializa el logger para escribir en archivo y consola
    static void initLogger() {
        try {
            tee = new PrintStream(new FileOutputStream(LOG_FILE, false), true, StandardCharsets.UTF_8);
        } catch (Exception e) {
            tee = null;
        }
        log("[" + TS.format(new Date()) + "] Logger abierto");
    }

    // cierra el logger si estaba abierto
    static void closeLogger() {
        if (tee != null) tee.close();
    }

    // escribe en consola y log file
    static void log(String s) {
        String msg = TS.format(new Date()) + " | " + s;
        System.out.println(msg);
        if (tee != null) tee.println(msg);
    }

    // expande las piezas según la cantidad (qty)
    static List<PartInstance> expandInstances(List<PartSpec> specs) {
        List<PartInstance> out = new ArrayList<>();
        for (PartSpec s : specs) {
            int q = Math.max(1, s.qty);
            for (int i=0;i<q;i++) out.add(new PartInstance(s, i));
        }
        return out;
    }

    // calcula los ángulos permitidos por pieza
    static Map<PartSpec,int[]> computeAngleOptions(List<PartSpec> specs) {
        Map<PartSpec,int[]> m = new HashMap<>();
        for (PartSpec s : specs) {
            int step = Math.max(1, s.rotationStepDeg);
            ArrayList<Integer> a = new ArrayList<>();
            for (int ang = 0; ang < 360; ang += step) a.add(ang);
            int[] arr = new int[a.size()];
            for (int i=0;i<a.size();i++) arr[i] = a.get(i);
            m.put(s, arr);
        }
        return m;
    }

    // colisión estricta: permite tocar bordes pero no solaparse
    static boolean collideStrict(Geometry a, Geometry b) {
        return a.intersects(b) && !a.touches(b);
    }

    // intenta mover una geometría paso a paso en dirección (stepx, stepy)
    static double[] tryNudge(Polygon sheetInset, Geometry gClear, List<Geometry> placedClear,
                             double stepx, double stepy, int steps)
    {
        double dx=0, dy=0;
        for (int i=0;i<steps;i++) {
            Geometry t = GeomUtils.translate(gClear, stepx, stepy);
            if (!sheetInset.covers(t)) break; // covers: se permite borde
            boolean inter=false;
            for (Geometry pg : placedClear) { 
                if (collideStrict(pg, t)) { inter=true; break; } 
            }
            if (inter) break;
            gClear = t; dx += stepx; dy += stepy;
        }
        // luego prueba en eje vertical (usa stepx como magnitud vertical según diseño original)
        for (int i=0;i<steps;i++) {
            Geometry t = GeomUtils.translate(gClear, 0, stepx);
            if (!sheetInset.covers(t)) break;
            boolean inter=false;
            for (Geometry pg : placedClear) { 
                if (collideStrict(pg, t)) { inter=true; break; } 
            }
            if (inter) break;
            gClear = t; dy += stepx;
        }
        return new double[]{dx, dy};
    }

    // calcula el score local de una pieza candidata
    static double localScore(DecodeResult dr, Geometry cc) {
        Envelope ce = cc.getEnvelopeInternal();
        double cy = ce.getMinY();
        double cx = ce.getMinX();

        double nGlobMinX = Math.min(dr.globMinX, ce.getMinX());
        double nGlobMinY = Math.min(dr.globMinY, ce.getMinY());
        double nGlobMaxX = Math.max(dr.globMaxX, ce.getMaxX());
        double nGlobMaxY = Math.max(dr.globMaxY, ce.getMaxY());
        double bboxGrow = (nGlobMaxX - nGlobMinX)*(nGlobMaxY - nGlobMinY)
                        - Math.max(0, (dr.globMaxX - dr.globMinX)*(dr.globMaxY - dr.globMinY));

        return W_BEST_Y*cy + W_BEST_X*cx + W_BBOX_INFLATE*bboxGrow;
    }

    // score usado en búsqueda local (solo bbox)
    static double localScoreBBoxOnly(Polygon sheet, DecodeResult dr, Geometry cc) {
        Envelope ce = cc.getEnvelopeInternal();
        double nGlobMinX = Math.min(dr.globMinX, ce.getMinX());
        double nGlobMinY = Math.min(dr.globMinY, ce.getMinY());
        double nGlobMaxX = Math.max(dr.globMaxX, ce.getMaxX());
        double nGlobMaxY = Math.max(dr.globMaxY, ce.getMaxY());
        return (nGlobMaxY - nGlobMinY) * (nGlobMaxX - nGlobMinX);
    }

    // inserta un candidato en la RCL manteniendo orden por score
    static void insertIntoRCL(List<Candidate> rcl, Candidate c, int k) {
        int i=0; 
        while (i<rcl.size() && rcl.get(i).scoreVal <= c.scoreVal) i++;
        rcl.add(i, c);
        if (rcl.size() > k) rcl.remove(rcl.size()-1);
    }
    // genera las anclas iniciales (esquinas + malla sobre la plancha)
    static ArrayList<Anchor> initialAnchors(double minX, double minY, double maxX, double maxY){
        ArrayList<Anchor> a = new ArrayList<>();
        // esquinas
        a.add(new Anchor(minX, minY));
        a.add(new Anchor(maxX, minY));
        a.add(new Anchor(minX, maxY));
        a.add(new Anchor(maxX, maxY));
        // malla repartida en toda la plancha
        for (int i=0; i<GRID_ANCHORS_X; i++) {
            for (int j=0; j<GRID_ANCHORS_Y; j++) {
                double x = minX + (i+0.5)*(maxX-minX)/GRID_ANCHORS_X;
                double y = minY + (j+0.5)*(maxY-minY)/GRID_ANCHORS_Y;
                a.add(new Anchor(x, y));
            }
        }
        return a;
    }

    // añade anclas locales alrededor de una geometría ya colocada
    static void addLocalAnchors(List<Anchor> anchors, Envelope ce, double minX, double minY, double maxX, double maxY){
        double rx = ce.getMaxX(), lx = ce.getMinX(), by = ce.getMinY(), ty = ce.getMaxY();
        double[] xs = {rx, lx};
        double[] ys = {by, ty};
        for (double x : xs) {
            anchors.add(new Anchor(x, by));
            anchors.add(new Anchor(x, ty));
        }
        for (double y : ys) {
            anchors.add(new Anchor(lx, y));
            anchors.add(new Anchor(rx, y));
        }
        anchors.removeIf(an -> an.x < minX-EPS_EDGE || an.x > maxX+EPS_EDGE || an.y < minY-EPS_EDGE || an.y > maxY+EPS_EDGE);
    }

    // selecciona aleatoriamente un subconjunto de anclas
    static List<Anchor> sampleAnchors(List<Anchor> anchors, Random rng, int max){
        if (anchors.size() <= max) return new ArrayList<>(anchors);
        ArrayList<Anchor> out = new ArrayList<>(max);
        ArrayList<Integer> idx = new ArrayList<>(anchors.size());
        for (int i=0;i<anchors.size();i++) idx.add(i);
        Collections.shuffle(idx, rng);
        for (int i=0;i<max;i++) out.add(anchors.get(idx.get(i)));
        return out;
    }

    // busca el PartSpec correspondiente a un id
    static PartSpec findSpec(Map<PartSpec,int[]> angleOptions, String id){
        for (PartSpec s : angleOptions.keySet()) if (s.id.equals(id)) return s;
        return null;
    }

    // amplía los ángulos con tweaks si el paso de rotación es pequeño
    static int[] maybeAugmentAngles(int[] base, PartSpec sp){
        if (sp.rotationStepDeg <= 5) {
            HashSet<Integer> set = new HashSet<>();
            for (int b : base) {
                set.add(norm(b));
                set.add(norm(b + TWEAK_ANGLE_STEP));
                set.add(norm(b - TWEAK_ANGLE_STEP));
            }
            int[] arr = new int[set.size()];
            int i=0; for (int v : set) arr[i++] = v;
            Arrays.sort(arr);
            return arr;
        }
        return base;
    }

    // normaliza ángulo a [0,360)
    static int norm(int ang) { 
        int a = ang % 360; 
        if (a < 0) a += 360; 
        return a; 
    }

    // calcula el área total ocupada por las piezas colocadas
    static double totalArea(List<PlacedPart> placed){
        double s=0; 
        for (PlacedPart p : placed) s += p.polyPlaced.getArea(); 
        return s;
    }

    // devuelve la coordenada Y máxima alcanzada por las piezas
    static double usedYmax(List<PlacedPart> placed){
        double y=0; 
        for (PlacedPart p : placed) 
            y = Math.max(y, p.polyPlaced.getEnvelopeInternal().getMaxY()); 
        return y;
    }

    // shuffle por bloques para variar el orden de piezas
    static <T> void blockShuffle(List<T> list, Random rng, int blockSize){
        for (int i=0; i<list.size(); i+=blockSize) {
            int end = Math.min(list.size(), i+blockSize);
            Collections.shuffle(list.subList(i, end), rng);
        }
    }

    // calcula coeficiente de variación de un array
    static double cv(double[] a){
        if (a==null || a.length==0) return 1.0;
        double s=0; 
        for(double x:a) s+=x;
        double m=s/a.length; 
        if (m==0) return 1.0;
        double v=0; 
        for(double x:a) v+=(x-m)*(x-m);
        v/=a.length; 
        return Math.sqrt(v)/m;
    }

    // calcula mediana de un array
    static double median(double[] a){
        if (a==null || a.length==0) return 0.0;
        double[] cp=a.clone(); 
        Arrays.sort(cp); 
        return cp[cp.length/2];
    }
    // llena la plancha con piezas usando rejilla (sin huecos, modo GRID fast-path)
    static void gridFillRaster(Polygon sheetInset,
                               List<PlacedPart> placed, List<Geometry> placedClear,
                               List<PartInstance> instances, Map<PartSpec,int[]> angleOptions,
                               double minX,double minY,double maxX,double maxY,
                               double cellW,double cellH)
    {
        final double EPS = 1e-9;
        int nx = Math.max(1, (int)Math.floor(((maxX - minX) + EPS) / cellW));
        int ny = Math.max(1, (int)Math.floor(((maxY - minY) + EPS) / cellH));

        Map<String,Integer> remaining = new HashMap<>();
        for (PartInstance pi : instances) remaining.merge(pi.spec.id, 1, Integer::sum);

        for (int r=0; r<ny; r++){
            double anchorY = minY + r*cellH;
            for (int c=0; c<nx; c++){
                double anchorX = minX + c*cellW;

                boolean placedHere = false;
                for (PartInstance pi : instances){
                    Integer rem = remaining.get(pi.spec.id);
                    if (rem==null || rem<=0) continue;

                    int[] allowed = angleOptions.getOrDefault(pi.spec, new int[]{0});
                    for (int angDeg: allowed){
                        Geometry g0  = GeomUtils.orient(pi.spec.poly, angDeg, false);
                        Geometry gc0 = (pi.spec.clearance > 0) ? g0.buffer(pi.spec.clearance, 8) : g0;

                        Envelope pe0 = gc0.getEnvelopeInternal();
                        double dxa = anchorX - pe0.getMinX();
                        double dya = anchorY - pe0.getMinY();

                        Geometry ccTry = GeomUtils.translate(gc0, dxa, dya);
                        if (!sheetInset.covers(ccTry)) continue;   // permite apoyo en borde

                        boolean clash = false;
                        for (Geometry pg : placedClear) {
                            if (collideStrict(pg, ccTry)) { clash = true; break; }
                        }
                        if (clash) continue;

                        Geometry cand = GeomUtils.translate(g0, dxa, dya);

                        PlacedPart pp = new PlacedPart();
                        pp.id = pi.spec.id; 
                        pp.polyPlaced = cand; 
                        pp.angleDeg = angDeg;
                        pp.mirrored = false; 
                        pp.dx = dxa; 
                        pp.dy = dya;

                        placed.add(pp);
                        placedClear.add(ccTry);
                        remaining.put(pi.spec.id, rem-1);
                        placedHere = true;
                        break;
                    }
                    if (placedHere) break;
                }
            }
        }
    }
}
