package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.prep.PreparedGeometry;
import org.locationtech.jts.geom.prep.PreparedGeometryFactory;
import org.locationtech.jts.index.quadtree.Quadtree;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/** MemeticPlacer mejorado: semillas BL-guided, insertion-moves, ruin&recreate,
 *  aprendizaje de ángulos, rates adaptativos por estancamiento,
 *  caches/PreparedPolygon/STRtree en decode, beam-lite con lookahead
 *  y fitness con penalización de perímetro expuesto. */
public class MemeticPlacer {

    // =================== Parámetros GA/Memético (tunable) ===================
    static final int    POP_SIZE          = 36;
    static final int    ELITE_SIZE        = 6;   // ↑ ligerísimo para reforzar buenas estructuras
    static final int    GEN_MAX           = 200;

    static final long   MAX_MILLIS_TOTAL  = 60_000;
    static final long   MAX_MILLIS_LS     = 20_000;

    static final double MUT_RATE_ORDER_BASE = 0.35;
    static final double MUT_RATE_ANGLE_BASE = 0.25;
    static final int    TOURNAMENT_K        = 4; // ↑ presión selectiva leve

    static final long   DEFAULT_SEED      = 42L;

    // =================== Parámetros de colocación ===================
    static final int GRID_ANCHORS_X      = 10;
    static final int GRID_ANCHORS_Y      = 6;
    static final int MAX_ANCHORS_PER_TRY = 220;
    static final int MAX_CHECKS_PER_PART = 900;   // ↓ menor para evitar demasiado trabajo por pieza
    static final double NUDGE_STEP       = 1.0;   // ↑ pasos más grandes, ya hay "gravity"
    static final int NUDGE_STEPS_MAX     = 6;

    static final double W_BEST_Y       = 1.0;
    static final double W_BEST_X       = 0.15;
    static final double W_BBOX_INFLATE = 0.01;

    static final double GRID_CV_MAX = 0.03; // ↑ activar grid-raster en más casos homogéneos
    static final double EPS_EDGE    = 1e-9;

    // ======= Accel: caches y estructuras rápidas =======
    static final boolean USE_PREPARED_SHEET = true;
    static final boolean USE_STR_TREE       = true;

    static class RotCacheEntry {
        final Geometry base;      // pieza rotada sin clearance
        final Geometry clear;     // base con clearance (buffer)
        final Envelope envClear;  // envelope del clear
        RotCacheEntry(Geometry base, Geometry clear){
            this.base = base; this.clear = clear; this.envClear = clear.getEnvelopeInternal();
        }
    }
    static final Map<String, RotCacheEntry> ROT_CACHE = new HashMap<>(4096);

    static String rotKey(PartSpec s, int angle){
        return s.id + "|" + angle; // ids de PartSpec deben ser únicos
    }
    static RotCacheEntry getRotCached(PartSpec s, int angle){
        String k = rotKey(s, angle);
        RotCacheEntry e = ROT_CACHE.get(k);
        if (e != null) return e;
        Geometry g0  = GeomUtils.orient(s.poly, angle, false);
        Geometry gc0 = (s.clearance > 0) ? g0.buffer(s.clearance, 8) : g0;
        e = new RotCacheEntry(g0, gc0);
        ROT_CACHE.put(k, e);
        return e;
    }

    // =================== Beam-lite / lookahead ===================
    static final int     BEAM_K               = 4;   // top-K anclas locales
    static final boolean USE_BEAM_LOOKAHEAD   = true;

    // =================== Métrica de fitness ===================
    static final double LAMBDA_PERIM = 0.002; 

    // =================== Logging ===================
    static final String LOG_FILE = "LogMemetic.log";
    static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");
    static final DecimalFormat DF = new DecimalFormat("#,##0.###");
    private static PrintStream tee;

    // ==== watchdogs/telemetría decode ====
    static final long  MAX_MILLIS_DECODE     = 15_000;
    static final int   LOG_DEC_POS_STEP      = 5;
    static final int   LOG_DEC_ANCHOR_STEP   = 80;
    static final int   GRAVITY_STEPS_MAX     = 600;

    // =================== Estructuras ===================
    static class PartInstance {
        final PartSpec spec; final int copyIdx;
        PartInstance(PartSpec s, int k) { spec = s; copyIdx = k; }
        String key() { return spec.id + "#" + copyIdx; }
    }
    static class Anchor { final double x, y; Anchor(double x, double y){ this.x=x; this.y=y; } }

    static class Candidate {
        Geometry clearGeom;
        Geometry placedGeom;
        double dx, dy;
        int angleDeg;
        double scoreVal;
    }

    static class DecodeResult {
        final List<PlacedPart> placed = new ArrayList<>();
        final List<Geometry>   placedClear = new ArrayList<>();
        double globMinX, globMinY, globMaxX, globMaxY;
    }

    static class Individual {
        int[] order;
        int[] angles;
        List<PlacedPart> layout;
        double area;     // área verdadera colocada (para reportes)
        double fitness;  // área - lambda * perímetro (para comparación)
        double yMax;
        long   tieNoise;
        long   fingerprint;
    }

    // =================== API pública ===================
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        long seed = randomSeed();
        return place(sheetInset, specs, seed);
    }

    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs, long seed) {
        initLogger();
        final long t0 = System.currentTimeMillis();

        List<PartInstance> instances = expandInstances(specs);
        final int N = instances.size();
        if (N == 0) { log("No hay instancias"); closeLogger(); return Collections.emptyList(); }

        Map<PartSpec,int[]> angleOptionsBySpec = computeAngleOptions(specs);
        int[][] angleDomain = new int[N][];
        for (int i=0;i<N;i++) {
            int[] dom = angleOptionsBySpec.get(instances.get(i).spec);
            angleDomain[i] = (dom == null || dom.length == 0) ? new int[]{0} : dom;
        }

        Random rng = new Random(seed);

        // cache cromosomas ya evaluados
        HashSet<Long> seen = new HashSet<>(POP_SIZE*4);

        List<Individual> pop = new ArrayList<>(POP_SIZE);

        // Semillas fuertes
        int greedyVariants = Math.min(ELITE_SIZE, POP_SIZE/3);
        for (int g=0; g<greedyVariants; g++) {
            Individual s = greedySeed(N, instances, angleDomain, rng);
            evaluate(s, sheetInset, instances, angleDomain, rng);
            registerFingerprint(s);
            if (seen.add(s.fingerprint)) pop.add(s);
        }

        // Semilla guiada por BL
        if (pop.size() < POP_SIZE) {
            Individual bl = blGuidedSeed(N, instances, angleDomain, sheetInset, rng);
            evaluate(bl, sheetInset, instances, angleDomain, rng);
            registerFingerprint(bl);
            if (seen.add(bl.fingerprint)) pop.add(bl);
        }

        // Relleno con aleatorios (evitando duplicados)
        while (pop.size() < POP_SIZE) {
            Individual ind = randomIndividual(N, angleDomain, rng);
            evaluate(ind, sheetInset, instances, angleDomain, rng);
            registerFingerprint(ind);
            if (seen.add(ind.fingerprint)) pop.add(ind);
        }

        Individual best = bestOf(pop);
        log("INI | seed=" + seed + " | N=" + N + " | pop=" + POP_SIZE +
            " | area=" + DF.format(best.area) + " | yMax=" + DF.format(best.yMax));

        // control adaptativo por estancamiento
        int stagnation = 0;
        double mutOrder = MUT_RATE_ORDER_BASE;
        double mutAngle = MUT_RATE_ANGLE_BASE;

        // histograma de ángulos (angle-learning)
        Map<Integer, Integer> angleHist = new HashMap<>();

        for (int gen = 0; gen < GEN_MAX; gen++) {
            if (System.currentTimeMillis() - t0 > MAX_MILLIS_TOTAL) {
                log("Tiempo total agotado en gen=" + gen); break;
            }

            pop.sort(Comparator.<Individual>comparingDouble((Individual a)->-a.fitness)
                               .thenComparingDouble(a -> a.yMax)
                               .thenComparingLong(a -> a.tieNoise));

            List<Individual> next = new ArrayList<>(POP_SIZE);
            for (int e=0; e<ELITE_SIZE && e<pop.size(); e++) next.add(cloneShallow(pop.get(e)));

            // refresco histograma de ángulos a partir de élites
            angleHist.clear();
            for (int e=0; e<Math.min(ELITE_SIZE, pop.size()); e++) {
                for (int a : pop.get(e).angles) angleHist.merge(a, 1, Integer::sum);
            }

            // construir hijos
            while (next.size() < POP_SIZE) {
                Individual p1 = tournament(pop, rng, TOURNAMENT_K);
                Individual p2 = tournament(pop, rng, TOURNAMENT_K);
                Individual child = crossoverOX_withAngleBlend(p1, p2, angleDomain, rng);

                adaptiveMutate(child, angleDomain, rng, mutOrder, mutAngle, angleHist);

                evaluate(child, sheetInset, instances, angleDomain, rng);
                registerFingerprint(child);
                if (seen.add(child.fingerprint)) next.add(child);

                // cada 4 hijos, añade un path-relinking (best -> elite aleatorio)
                if (next.size() < POP_SIZE && next.size() % 4 == 0) {
                    Individual anchor = next.get(0); // mejor actual entre los que ya están en next
                    Individual elite2 = pop.get(rng.nextInt(Math.min(ELITE_SIZE, pop.size())));
                    Individual pr = pathRelink(anchor, elite2, angleDomain, rng);
                    evaluate(pr, sheetInset, instances, angleDomain, rng);
                    registerFingerprint(pr);
                    if (seen.add(pr.fingerprint)) next.add(pr);
                }
            }

            pop = next;

            // ===================== Búsqueda local =====================
            long lsBudget = Math.max(0, MAX_MILLIS_LS - (System.currentTimeMillis() - t0) / 6);
            if (lsBudget > 0) {
                int elites = Math.min(ELITE_SIZE, pop.size());
                long perElite = Math.max(5L, lsBudget / elites);

                for (int e=0; e<elites; e++) {
                    localSearchStrong(pop.get(e), sheetInset, instances, angleDomain, rng, perElite);
                }
            }

            Individual curBest = bestOf(pop);
            if (isBetter(curBest, best)) {
                best = cloneDeep(curBest);
                stagnation = 0;
                mutOrder = Math.max(0.15, mutOrder * 0.92);
                mutAngle = Math.max(0.12, mutAngle * 0.92);
            } else {
                stagnation++;
                if (stagnation % 8 == 0) {
                    mutOrder = Math.min(0.65, mutOrder * 1.15);
                    mutAngle = Math.min(0.55, mutAngle * 1.15);
                }
            }

            if (gen>0 && gen % 25 == 0) {
                pop.sort(Comparator.<Individual>comparingDouble(a -> a.fitness)
                                   .thenComparingDouble(a -> -a.yMax));
                int pert = Math.max(2, POP_SIZE/6);
                for (int i=0; i<pert; i++) {
                    Individual w = pop.get(i);
                    weakPerturb(w, angleDomain, rng);
                    evaluate(w, sheetInset, instances, angleDomain, rng);
                }
            }

            log(String.format("GEN %d | best area=%s yMax=%s | placed=%d | fitness=%.3f | mutOrder=%.2f mutAngle=%.2f",
                    gen, DF.format(best.area), DF.format(best.yMax),
                    (best.layout==null?0:best.layout.size()), best.fitness, mutOrder, mutAngle));
        }

        int totalPiezas        = instances.size();
        int colocadas          = (best.layout == null ? 0 : best.layout.size());
        double areaPlancha     = sheetInset.getArea();
        double areaColocada    = best.area; // área verdadera (no penalizada)
        double areaDesperd     = Math.max(0.0, areaPlancha - areaColocada);
        double aprovechamiento = areaPlancha > 0 ? (100.0 * areaColocada / areaPlancha) : 0.0;
        long durTotalMs        = System.currentTimeMillis() - t0;

        log("------------------------------");
        log("RESUMEN FINAL (MEMETIC - single run)");
        log("Seed               : " + seed);
        log("Piezas a posicionar: " + totalPiezas);
        log("Piezas colocadas   : " + colocadas);
        log("% aprovechamiento  : " + DF.format(aprovechamiento) + " %");
        log("Área colocada      : " + DF.format(areaColocada));
        log("Área plancha       : " + DF.format(areaPlancha));
        log("Área desperdiciada : " + DF.format(areaDesperd));
        log("Ymax               : " + DF.format(best.yMax));
        log("Tiempo total       : " + durTotalMs + " ms (" + DF.format(durTotalMs/1000.0) + " s)");
        log("------------------------------");
        log("BEST | colocadas=" + colocadas + " | area=" + DF.format(areaColocada) + " | yMax=" + DF.format(best.yMax));

        closeLogger();
        return best.layout == null ? Collections.emptyList() : best.layout;
    }

    // =================== Decode (con caches/Prepared/STRtree y beam-lite) ===================
    static DecodeResult decode(Polygon sheetInset,
                           List<PartInstance> orderInst,
                           int[] angles,
                           Random rng)
    {
        final long tStart = System.currentTimeMillis();
        DecodeResult res = new DecodeResult();
        Envelope env = sheetInset.getEnvelopeInternal();
        double minX = env.getMinX(), minY = env.getMinY(), maxX = env.getMaxX(), maxY = env.getMaxY();
        res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;

        log("decode | start | parts=" + orderInst.size());

        // Prepared geometry para covers/contains
        final PreparedGeometry prepSheet = USE_PREPARED_SHEET ? PreparedGeometryFactory.prepare(sheetInset) : null;
        // Índice espacial para colisiones
        final Quadtree qt = USE_STR_TREE ? new Quadtree() : null;

        boolean GRID_MODE = false;
        double cellW = 0, cellH = 0;
        {
            double rectArea = (maxX - minX) * (maxY - minY);
            boolean sheetRectangular = Math.abs(sheetInset.getArea() - rectArea) < 1e-6;
            int limit = Math.min(50, orderInst.size());
            if (sheetRectangular && limit > 0) {
                double[] ws = new double[limit], hs = new double[limit];
                for (int i=0;i<limit;i++){
                    PartSpec sp = orderInst.get(i).spec;
                    int angDeg = (angles.length>i)? angles[i] : 0;
                    RotCacheEntry rc = getRotCached(sp, angDeg);
                    Envelope e   = rc.envClear;
                    ws[i] = e.getWidth(); hs[i] = e.getHeight();
                }
                if (cv(ws) <= GRID_CV_MAX && cv(hs) <= GRID_CV_MAX) {
                    GRID_MODE = true; cellW = median(ws); cellH = median(hs);
                    log("decode | GRID_MODE=true | cellW=" + DF.format(cellW) + " cellH=" + DF.format(cellH));
                }
            }
        }

        if (GRID_MODE) {
            gridFillRaster(sheetInset, res.placed, res.placedClear, orderInst, angles,
                           minX, minY, maxX, maxY, cellW, cellH);
            // bbox global
            res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;
            for (Geometry cc : res.placedClear) {
                Envelope ce = cc.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());
                if (USE_STR_TREE) qt.insert(ce, cc);
            }
            log("decode | end (grid) | placed=" + res.placed.size() +
                " | elapsedMs=" + (System.currentTimeMillis()-tStart));
            return res;
        }

        ArrayList<Anchor> anchors = initialAnchors(minX, minY, maxX, maxY);

        for (int pos = 0; pos < orderInst.size(); pos++) {

            // --- Watchdog de tiempo total del decode ---
            long elapsed = System.currentTimeMillis() - tStart;
            if (elapsed > MAX_MILLIS_DECODE) {
                log("decode | TIMEOUT | pos=" + pos + "/" + orderInst.size() +
                    " | placed=" + res.placed.size() +
                    " | elapsedMs=" + elapsed + " | aborting decode (partial layout returned)");
                return res; // devolvemos lo que haya
            }

            if (pos % LOG_DEC_POS_STEP == 0) {
                log("decode | progress | pos=" + pos + "/" + orderInst.size() +
                    " | placed=" + res.placed.size() +
                    " | elapsedMs=" + elapsed);
            }

            PartInstance inst = orderInst.get(pos);
            int angDeg = norm(angles[pos]);

            RotCacheEntry rc = getRotCached(inst.spec, angDeg);
            Geometry gBase = rc.base;
            Geometry gClr  = rc.clear;
            Envelope pe    = rc.envClear;

            List<Anchor> sampleAnchors = sampleAnchors(anchors, rng, MAX_ANCHORS_PER_TRY);
            // Beam-lite: mantenemos top-K candidatos
            PriorityQueue<Candidate> topK = new PriorityQueue<>(Comparator.comparingDouble(c -> -c.scoreVal)); // max-heap por score
            int checks=0;

            for (int idxA = 0; idxA < sampleAnchors.size(); idxA++) {
                Anchor a = sampleAnchors.get(idxA);

                if (idxA > 0 && idxA % LOG_DEC_ANCHOR_STEP == 0) {
                    long elA = System.currentTimeMillis() - tStart;
                    log("decode | pos=" + pos + " | anchorsTried=" + idxA +
                        "/" + sampleAnchors.size() + " | placed=" + res.placed.size() +
                        " | elapsedMs=" + elA);
                }

                if (++checks > MAX_CHECKS_PER_PART) break;

                double dxa = a.x - pe.getMinX();
                double dya = a.y - pe.getMinY();

                double txMin = pe.getMinX() + dxa, txMax = pe.getMaxX() + dxa;
                double tyMin = pe.getMinY() + dya, tyMax = pe.getMaxY() + dya;
                if (txMin < minX-EPS_EDGE || txMax > maxX+EPS_EDGE || tyMin < minY-EPS_EDGE || tyMax > maxY+EPS_EDGE) continue;

                Geometry cc = GeomUtils.translate(gClr, dxa, dya);
                boolean okInSheet = USE_PREPARED_SHEET ? prepSheet.covers(cc) : sheetInset.covers(cc);
                if (!okInSheet) continue;

                boolean inter = false;
                if (USE_STR_TREE && !res.placedClear.isEmpty()) {
                    @SuppressWarnings("unchecked")
                    List<Geometry> cand = qt.query(cc.getEnvelopeInternal());
                    for (Geometry pg : cand) { if (collideStrict(pg, cc)) { inter = true; break; } }
                } else {
                    for (Geometry pg : res.placedClear) { if (collideStrict(pg, cc)) { inter = true; break; } }
                }
                if (inter) continue;

                // nudge primario
                double[] nudged = tryNudge(sheetInset, cc, res.placedClear, -NUDGE_STEP, -NUDGE_STEP, NUDGE_STEPS_MAX);
                cc = GeomUtils.translate(cc, nudged[0], nudged[1]);
                Geometry placed = GeomUtils.translate(gBase, dxa+nudged[0], dya+nudged[1]);

                // gravity con tope de pasos
                double[] grav = gravityCompactCapped(sheetInset, cc, res.placedClear, GRAVITY_STEPS_MAX);
                cc = GeomUtils.translate(cc, grav[0], grav[1]);
                placed = GeomUtils.translate(placed, grav[0], grav[1]);

                double val = localScore(res, cc);

                Candidate candC = new Candidate();
                candC.clearGeom=cc; candC.placedGeom=placed;
                candC.dx = dxa+nudged[0]+grav[0]; candC.dy = dya+nudged[1]+grav[1];
                candC.angleDeg = angDeg; candC.scoreVal = val;

                topK.offer(candC);
                if (topK.size() > BEAM_K) topK.poll(); // mantenemos K mejores (heap elimina el peor)
            }

            Candidate chosen = null;
            if (!topK.isEmpty()) {
                if (USE_BEAM_LOOKAHEAD && pos+1 < orderInst.size()) {
                    chosen = pickWithLookahead(topK, orderInst, angles, pos+1, res, rng);
                } else {
                    chosen = topK.stream().min(Comparator.comparingDouble(c -> c.scoreVal)).get();
                }
            }

            if (chosen != null) {
                PlacedPart pp = new PlacedPart();
                pp.id = orderInst.get(pos).spec.id;
                pp.polyPlaced = chosen.placedGeom;
                pp.angleDeg = chosen.angleDeg;
                pp.mirrored = false;
                pp.dx = chosen.dx; pp.dy = chosen.dy;

                res.placed.add(pp);
                res.placedClear.add(chosen.clearGeom);
                if (USE_STR_TREE) qt.insert(chosen.clearGeom.getEnvelopeInternal(), chosen.clearGeom);

                Envelope ce = chosen.clearGeom.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());

                addLocalAnchors(anchors, ce, minX, minY, maxX, maxY);
            }
        }

        log("decode | end | placed=" + res.placed.size() +
            " | elapsedMs=" + (System.currentTimeMillis()-tStart));
        return res;
    }

    // ==== lookahead de profundidad 1 (barato, sobre envelopes) ====
    static Candidate pickWithLookahead(
            PriorityQueue<Candidate> topK,
            List<PartInstance> orderInst, int[] angles, int nextPos,
            DecodeResult partial, Random rng)
    {
        double gMinX = partial.globMinX, gMinY = partial.globMinY, gMaxX = partial.globMaxX, gMaxY = partial.globMaxY;

        Candidate bestC = null;
        double bestScore = Double.POSITIVE_INFINITY;

        List<Candidate> cands = new ArrayList<>(topK);
        for (Candidate c : cands) {
            // aplicar efectos bbox del actual
            Envelope ce = c.clearGeom.getEnvelopeInternal();
            double nMinX = Math.min(gMinX, ce.getMinX());
            double nMinY = Math.min(gMinY, ce.getMinY());
            double nMaxX = Math.max(gMaxX, ce.getMaxX());
            double nMaxY = Math.max(gMaxY, ce.getMaxY());

            double nextPenalty = 0.0;
            if (nextPos < orderInst.size()) {
                PartInstance nxt = orderInst.get(nextPos);
                int angDeg = norm(angles[nextPos]);

                RotCacheEntry rc = getRotCached(nxt.spec, angDeg);
                Envelope npe = rc.envClear;

                double[][] estAnchors = new double[][]{
                        {nMinX, nMinY}, {nMaxX, nMinY}, {nMinX, nMaxY}, {nMaxX, nMaxY}
                };
                double bestLocal = Double.POSITIVE_INFINITY;
                for (double[] a : estAnchors) {
                    double dxa = a[0] - npe.getMinX();
                    double dya = a[1] - npe.getMinY();
                    Envelope tryEnv = new Envelope(npe);
                    tryEnv.translate(dxa, dya);

                    double nnMinX = Math.min(nMinX, tryEnv.getMinX());
                    double nnMinY = Math.min(nMinY, tryEnv.getMinY());
                    double nnMaxX = Math.max(nMaxX, tryEnv.getMaxX());
                    double nnMaxY = Math.max(nMaxY, tryEnv.getMaxY());
                    double grow = (nnMaxX-nnMinX)*(nnMaxY-nnMinY) - (nMaxX-nMinX)*(nMaxY-nMinY);

                    double candScore = W_BEST_Y*tryEnv.getMinY() + W_BEST_X*tryEnv.getMinX() + W_BBOX_INFLATE*grow;
                    if (candScore < bestLocal) bestLocal = candScore;
                }
                nextPenalty = bestLocal;
            }

            double combined = c.scoreVal + 0.25 * nextPenalty; // peso suave
            if (combined < bestScore) { bestScore = combined; bestC = c; }
        }
        return bestC;
    }

    // ==== variante capada de gravity ====
    static double[] gravityCompactCapped(Polygon sheet, Geometry gClear, List<Geometry> placedClear, int maxSteps) {
        double dx=0, dy=0;
        int steps = 0;

        // bajar
        while (steps < maxSteps) {
            Geometry t = GeomUtils.translate(gClear, 0, -NUDGE_STEP);
            if (!sheet.covers(t)) break;
            boolean inter=false; for (Geometry pg:placedClear){ if (collideStrict(pg,t)){ inter=true; break; } }
            if (inter) break;
            gClear = t; dy -= NUDGE_STEP; steps++;
        }
        // izquierda
        while (steps < maxSteps) {
            Geometry t = GeomUtils.translate(gClear, -NUDGE_STEP, 0);
            if (!sheet.covers(t)) break;
            boolean inter=false; for (Geometry pg:placedClear){ if (collideStrict(pg,t)){ inter=true; break; } }
            if (inter) break;
            gClear = t; dx -= NUDGE_STEP; steps++;
        }
        return new double[]{dx, dy};
    }

    // =================== Búsqueda local fuerte ===================
    static void localSearchStrong(Individual ind,
                                  Polygon sheetInset,
                                  List<PartInstance> instances,
                                  int[][] angleDomain,
                                  Random rng,
                                  long millisBudget)
    {
        long t0 = System.currentTimeMillis();
        if (ind.layout == null || ind.layout.size() <= 1) return;

        // (1) Tweaks de ángulos
        int tries = Math.min(20, ind.layout.size());
        for (int t = 0; t < tries; t++) {
            if (System.currentTimeMillis() - t0 > millisBudget) return;

            int idxGene = rng.nextInt(ind.order.length);
            int[] dom = angleDomain[ind.order[idxGene]];
            if (dom != null && dom.length >= 2) {
                int cur = snapToDomain(ind.angles[idxGene], dom);
                int pos = Arrays.binarySearch(dom, cur);
                if (pos < 0) pos = -pos - 1;
                int step = rng.nextBoolean() ? 1 : -1;
                ind.angles[idxGene] = dom[(pos + step + dom.length) % dom.length];

                if (dom.length >= 36 && rng.nextDouble() < 0.6) {
                    int extra = 1 + rng.nextInt(Math.min(5, dom.length - 1));
                    int dir = rng.nextBoolean() ? 1 : -1;
                    int p2 = Arrays.binarySearch(dom, ind.angles[idxGene]);
                    if (p2 < 0) p2 = -p2 - 1;
                    ind.angles[idxGene] = dom[(p2 + dir * extra + dom.length) % dom.length];
                }
            }
            evaluate(ind, sheetInset, instances, angleDomain, rng);
        }

        // (2) Insertion-moves
        int moves = Math.min(30, ind.order.length);
        for (int m=0; m<moves; m++) {
            if (System.currentTimeMillis() - t0 > millisBudget) return;
            int i = rng.nextInt(ind.order.length);
            int j = rng.nextInt(ind.order.length);
            if (i == j) continue;
            insertion(ind.order, i, j);
            rotateAngleWithGene(ind, i, j);
            evaluate(ind, sheetInset, instances, angleDomain, rng);
        }

        // (3) Ruin & Recreate ligero
        if (System.currentTimeMillis() - t0 <= millisBudget) {
            int len = Math.max(2, ind.order.length/10);
            int s = rng.nextInt(Math.max(1, ind.order.length - len));
            ruinAndRecreate(ind, s, len, angleDomain, rng);
            evaluate(ind, sheetInset, instances, angleDomain, rng);
        }
    }

    // =================== Mutación adaptativa con “angle-learning” ===================
    static void adaptiveMutate(Individual ch, int[][] angleDomain, Random rng,
                               double mutOrder, double mutAngle,
                               Map<Integer,Integer> angleHist)
    {
        // orden
        if (rng.nextDouble() < mutOrder) {
            int i = rng.nextInt(ch.order.length);
            int j = rng.nextInt(ch.order.length);
            swap(ch.order, i, j);
            int tmp = ch.angles[i]; ch.angles[i] = ch.angles[j]; ch.angles[j] = tmp;
        }
        if (rng.nextDouble() < mutOrder*0.75) {
            int i = rng.nextInt(ch.order.length-1);
            reverse(ch.order, i, Math.min(ch.order.length-1, i+1+rng.nextInt(Math.max(1, ch.order.length/6))));
        }
        if (rng.nextDouble() < mutOrder*0.35) {
            int i = rng.nextInt(ch.order.length);
            int j = rng.nextInt(ch.order.length);
            insertion(ch.order, i, j);
            rotateAngleWithGene(ch, i, j);
        }

        // ángulos con sesgo por histograma
        for (int i=0;i<ch.angles.length;i++) {
            if (rng.nextDouble() < mutAngle) {
                int gene = ch.order[i];
                int[] dom = angleDomain[gene];
                if (dom.length == 1) { ch.angles[i] = dom[0]; continue; }

                if (rng.nextDouble() < 0.6) {
                    int idx = Arrays.binarySearch(dom, snapToDomain(ch.angles[i], dom));
                    if (idx < 0) idx = -idx-1;
                    int step = rng.nextBoolean()?1:-1;
                    ch.angles[i] = dom[(idx + step + dom.length) % dom.length];
                } else {
                    ch.angles[i] = sampleAngleByHist(dom, angleHist, rng);
                }
            }
        }
    }

    static int sampleAngleByHist(int[] dom, Map<Integer,Integer> hist, Random rng) {
        int tot=0; for (int a:dom) tot += Math.max(1, hist.getOrDefault(a, 1));
        int r = rng.nextInt(Math.max(1, tot));
        int acc=0;
        for (int a:dom) {
            acc += Math.max(1, hist.getOrDefault(a, 1));
            if (r < acc) return a;
        }
        return dom[rng.nextInt(dom.length)];
    }

    // =================== Evaluación / Comparadores ===================
    static void evaluate(Individual ind,
                         Polygon sheetInset,
                         List<PartInstance> instances,
                         int[][] angleDomain,
                         Random rng)
    {
        ArrayList<PartInstance> ord = new ArrayList<>(ind.order.length);
        for (int i=0;i<ind.order.length;i++) ord.add(instances.get(ind.order[i]));

        int[] angs = new int[ind.order.length];
        for (int i=0;i<angs.length;i++) angs[i] = ind.angles[i];

        DecodeResult dr = decode(sheetInset, ord, angs, rng);
        ind.layout = dr.placed;

        ind.area   = totalArea(dr.placed); // área verdadera reportable
        double perim = 0.0;
        for (PlacedPart p : dr.placed) perim += p.polyPlaced.getLength();
        ind.fitness = ind.area - LAMBDA_PERIM * perim;

        ind.yMax   = usedYmax(dr.placed);
    }

    static boolean isBetter(Individual a, Individual b){
        if ((a.fitness > b.fitness)) return true;
        if ((a.fitness < b.fitness)) return false;
        if (Math.abs(a.fitness-b.fitness) < 1e-9) {
            if (a.yMax < b.yMax) return true;
            if (a.yMax > b.yMax) return false;
            return Long.compare(a.tieNoise, b.tieNoise) < 0;
        }
        return false;
    }

    static Individual bestOf(List<Individual> pop){
        Individual best = pop.get(0);
        for (int i=1;i<pop.size();i++) if (isBetter(pop.get(i), best)) best = pop.get(i);
        return best;
    }

    // =================== Path-Relinking (best -> elite) ===================
    static Individual pathRelink(Individual a, Individual b, int[][] angleDomain, Random rng){
        Individual cur = cloneDeep(a);
        int n = cur.order.length;
        int[] posB = new int[n];
        for (int i=0;i<n;i++) posB[b.order[i]] = i;

        int steps = Math.min(n/4, 20);
        for (int s=0; s<steps; s++) {
            int pick = -1;
            for (int i=0;i<n;i++){ if (posB[cur.order[i]] != i){ pick=i; break; } }
            if (pick<0) break;
            int to = posB[cur.order[pick]];
            insertion(cur.order, pick, to);
            rotateAngleWithGene(cur, pick, to);

            int g = cur.order[to];
            if (rng.nextDouble()<0.35) {
                int ang = b.angles[posB[g]];
                cur.angles[to] = snapToDomain(ang, angleDomain[g]);
            }
        }
        return cur;
    }

    // =================== Semillas ===================
    static Individual randomIndividual(int N, int[][] angleDomain, Random rng){
        Individual ind = new Individual();
        ind.order = new int[N];
        for (int i=0;i<N;i++) ind.order[i] = i;
        shuffle(ind.order, rng);
        ind.angles = new int[N];
        for (int i=0;i<N;i++) {
            int[] dom = angleDomain[i];
            ind.angles[i] = dom[rng.nextInt(dom.length)];
        }
        ind.tieNoise = rng.nextLong();
        return ind;
    }

    static Individual greedySeed(int N, List<PartInstance> instances, int[][] angleDomain, Random rng){
        Integer[] idx = new Integer[N];
        for (int i=0;i<N;i++) idx[i]=i;

        // Precomputar áreas y un jitter fijo por elemento (consistente durante TODO el sort)
        double[] area = new double[N];
        double[] jitter = new double[N];
        for (int i=0;i<N;i++){
            double a = instances.get(i).spec.poly.getArea();
            area[i] = a;
            // jitter muy pequeño, escalado suavemente para no afectar el orden salvo empates
            jitter[i] = (rng.nextDouble() - 0.5) * 1e-9 * Math.max(1.0, a);
        }

        // Orden descendente por (area + jitter) — sin RNG dentro del comparador
        Arrays.sort(idx, (a,b) -> {
            double ka = area[a] + jitter[a];
            double kb = area[b] + jitter[b];
            int cmp = Double.compare(kb, ka); // desc
            if (cmp != 0) return cmp;
            // Tiebreaker determinista por índice para garantizar contrato total
            return Integer.compare(a, b);
        });

        Individual ind = new Individual();
        ind.order = new int[N]; 
        for (int i=0;i<N;i++) ind.order[i] = idx[i];

        ind.angles = new int[N];
        for (int i=0;i<N;i++) {
            int gene = ind.order[i];
            int[] dom = angleDomain[gene];
            ind.angles[i] = dom[rng.nextInt(dom.length)];
        }

        if (N >= 4 && rng.nextDouble()<0.5) {
            int start = rng.nextInt(N-2);
            reverse(ind.order, start, Math.min(N-1, start+2+rng.nextInt(3)));
        }
        ind.tieNoise = rng.nextLong();
        return ind;
    }


    /** Semilla “BL-guided”: ordena por altura y luego anchura del bbox con mejor ángulo discreto */
    static Individual blGuidedSeed(int N, List<PartInstance> instances, int[][] angleDomain, Polygon sheet, Random rng) {
        Integer[] idx = new Integer[N];
        for (int i=0;i<N;i++) idx[i]=i;

        double[] height = new double[N];
        double[] width  = new double[N];
        int[] bestAng   = new int[N];

        for (int i=0;i<N;i++) {
            int[] dom = angleDomain[i];
            double bestH = Double.POSITIVE_INFINITY, bestW = Double.POSITIVE_INFINITY;
            int    bestA = dom[0];
            for (int a : dom) {
                Geometry g = GeomUtils.orient(instances.get(i).spec.poly, a, false);
                if (instances.get(i).spec.clearance > 0) g = g.buffer(instances.get(i).spec.clearance, 8);
                Envelope e = g.getEnvelopeInternal();
                double h = e.getHeight(), w = e.getWidth();
                if (h < bestH || (Math.abs(h-bestH)<1e-9 && w < bestW)) {
                    bestH=h; bestW=w; bestA=a;
                }
            }
            height[i]=bestH; width[i]=bestW; bestAng[i]=bestA;
        }

        Arrays.sort(idx, (a,b)->{
            int c = Double.compare(height[b], height[a]); // alto desc
            if (c!=0) return c;
            return Double.compare(width[b], width[a]);    // ancho desc
        });

        Individual ind = new Individual();
        ind.order = new int[N];
        ind.angles = new int[N];
        for (int i=0;i<N;i++){
            ind.order[i]=idx[i];
            ind.angles[i]=bestAng[idx[i]];
        }
        ind.tieNoise = rng.nextLong();
        return ind;
    }

    // =================== Selección / Cruce ===================
    static Individual tournament(List<Individual> pop, Random rng, int k){
        Individual best = null;
        for (int i=0;i<k;i++) {
            Individual c = pop.get(rng.nextInt(pop.size()));
            if (best==null || isBetter(c, best)) best = c;
        }
        return best;
    }

    static Individual crossoverOX_withAngleBlend(Individual p1, Individual p2, int[][] angleDomain, Random rng){
        int N = p1.order.length;
        Individual ch = new Individual();
        ch.order = ox(p1.order, p2.order, rng);
        ch.angles = new int[N];

        boolean[] fromP1 = originMaskOX(p1.order, p2.order, ch.order);
        for (int i=0;i<N;i++) {
            int gene = ch.order[i];
            int a1 = p1.angles[indexOf(p1.order, gene)];
            int a2 = p2.angles[indexOf(p2.order, gene)];
            ch.angles[i] = fromP1[i] ? a1 : (rng.nextBoolean()?a1:a2);
            int[] dom = angleDomain[gene];
            ch.angles[i] = snapToDomain(ch.angles[i], dom);
        }
        ch.tieNoise = rng.nextLong();
        return ch;
    }

    // =================== Operadores auxiliares (orden/ángulos) ===================
    static void insertion(int[] a, int i, int j){
        if (i==j) return;
        int tmp = a[i];
        if (i<j) { System.arraycopy(a, i+1, a, i, j-i); a[j]=tmp; }
        else     { System.arraycopy(a, j, a, j+1, i-j); a[j]=tmp; }
    }

    static void rotateAngleWithGene(Individual ind, int from, int to){
        if (from==to) return;
        int tmp = ind.angles[from];
        if (from<to) { System.arraycopy(ind.angles, from+1, ind.angles, from, to-from); ind.angles[to]=tmp; }
        else         { System.arraycopy(ind.angles, to,   ind.angles, to+1, from-to);   ind.angles[to]=tmp; }
    }

    static void ruinAndRecreate(Individual ind, int start, int len, int[][] angleDomain, Random rng){
        int n = ind.order.length;
        int[] block = Arrays.copyOfRange(ind.order, start, start+len);
        int[] blockAngles = Arrays.copyOfRange(ind.angles, start, start+len);

        int[] rest = new int[n-len];
        int[] restAngles = new int[n-len];
        int p=0;
        for (int i=0;i<n;i++) if (i<start || i>=start+len) { rest[p]=ind.order[i]; restAngles[p]=ind.angles[i]; p++; }

        Integer[] orderB = new Integer[len];
        for (int i=0;i<len;i++) orderB[i]=i;
        Arrays.sort(orderB, (i,j)->Integer.compare(block[i], block[j])); // proxy simple

        List<Integer> newOrder = new ArrayList<>(n);
        List<Integer> newAngles= new ArrayList<>(n);
        for (int i=0;i<rest.length;i++){ newOrder.add(rest[i]); newAngles.add(restAngles[i]); }

        for (int id : orderB) {
            int pos = rng.nextInt(newOrder.size()+1);
            newOrder.add(pos, block[id]);
            newAngles.add(pos, blockAngles[id]);
        }

        for (int i=0;i<n;i++){ ind.order[i]=newOrder.get(i); ind.angles[i]=newAngles.get(i); }
    }

    static void weakPerturb(Individual w, int[][] angleDomain, Random rng){
        int len = Math.max(2, w.order.length/8);
        int s = rng.nextInt(Math.max(1, w.order.length-len));
        reverse(w.order, s, s+len-1);
        for (int k=0; k<w.angles.length; k++) {
            if (rng.nextDouble() < 0.15) {
                int[] dom = angleDomain[w.order[k]];
                w.angles[k] = dom[rng.nextInt(dom.length)];
            }
        }
    }

    // =================== Puntuación / Métricas ===================
    static boolean collideStrict(Geometry a, Geometry b) {
        return a.intersects(b) && !a.touches(b);
    }

    static double[] tryNudge(Polygon sheetInset, Geometry gClear, List<Geometry> placedClear,
                             double stepx, double stepy, int steps)
    {
        double dx=0, dy=0;
        for (int i=0;i<steps;i++) {
            Geometry t = GeomUtils.translate(gClear, 0, stepy);
            if (!sheetInset.covers(t)) break;
            boolean inter=false; for (Geometry pg:placedClear){ if (collideStrict(pg,t)){inter=true;break;} }
            if (inter) break; gClear=t; dy+=stepy;
        }
        for (int i=0;i<steps;i++) {
            Geometry t = GeomUtils.translate(gClear, stepx, 0);
            if (!sheetInset.covers(t)) break;
            boolean inter=false; for (Geometry pg:placedClear){ if (collideStrict(pg,t)){inter=true;break;} }
            if (inter) break; gClear=t; dx+=stepx;
        }
        for (int i=0;i<steps;i++) {
            Geometry t = GeomUtils.translate(gClear, stepx*0.5, stepy*0.5);
            if (!sheetInset.covers(t)) break;
            boolean inter=false; for (Geometry pg:placedClear){ if (collideStrict(pg,t)){inter=true;break;} }
            if (inter) break; gClear=t; dx+=stepx*0.5; dy+=stepy*0.5;
        }
        return new double[]{dx, dy};
    }

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

    static double totalArea(List<PlacedPart> placed){
        double s=0; for (PlacedPart p : placed) s += p.polyPlaced.getArea(); return s;
    }

    static double usedYmax(List<PlacedPart> placed){
        double y=0; for (PlacedPart p : placed) y = Math.max(y, p.polyPlaced.getEnvelopeInternal().getMaxY()); return y;
    }

    // =================== Anclas / rejilla ===================
    static ArrayList<Anchor> initialAnchors(double minX, double minY, double maxX, double maxY){
        ArrayList<Anchor> a = new ArrayList<>();
        a.add(new Anchor(minX, minY));
        a.add(new Anchor(maxX, minY));
        a.add(new Anchor(minX, maxY));
        a.add(new Anchor(maxX, maxY));
        for (int i=0; i<GRID_ANCHORS_X; i++) {
            for (int j=0; j<GRID_ANCHORS_Y; j++) {
                double x = minX + (i+0.5)*(maxX-minX)/GRID_ANCHORS_X;
                double y = minY + (j+0.5)*(maxY-minY)/GRID_ANCHORS_Y;
                a.add(new Anchor(x, y));
            }
        }
        return a;
    }

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

    static List<Anchor> sampleAnchors(List<Anchor> anchors, Random rng, int max){
        if (anchors.size() <= max) return new ArrayList<>(anchors);
        ArrayList<Anchor> out = new ArrayList<>(max);
        ArrayList<Integer> idx = new ArrayList<>(anchors.size());
        for (int i=0;i<anchors.size();i++) idx.add(i);
        Collections.shuffle(idx, rng);
        for (int i=0;i<max;i++) out.add(anchors.get(idx.get(i)));
        return out;
    }

    static void gridFillRaster(Polygon sheetInset,
                               List<PlacedPart> placed, List<Geometry> placedClear,
                               List<PartInstance> instances, int[] angles,
                               double minX,double minY,double maxX,double maxY,
                               double cellW,double cellH)
    {
        final double EPS = 1e-9;
        int nx = Math.max(1, (int)Math.floor(((maxX - minX) + EPS) / cellW));
        int ny = Math.max(1, (int)Math.floor(((maxY - minY) + EPS) / cellH));

        Map<String,Integer> remaining = new HashMap<>();
        for (PartInstance pi : instances) remaining.merge(pi.spec.id, 1, Integer::sum);

        List<PartInstance> ord = new ArrayList<>(instances);
        ord.sort((a,b)->Double.compare(b.spec.poly.getArea(), a.spec.poly.getArea()));

        for (int r=0; r<ny; r++){
            double anchorY = minY + r*cellH;
            for (int c=0; c<nx; c++){
                double anchorX = minX + c*cellW;

                for (PartInstance pi : ord){
                    Integer rem = remaining.get(pi.spec.id);
                    if (rem==null || rem<=0) continue;

                    int angDeg = 0;
                    int idx = instances.indexOf(pi);
                    if (idx >= 0 && idx < angles.length) angDeg = norm(angles[idx]);

                    RotCacheEntry rc = getRotCached(pi.spec, angDeg);
                    Envelope pe0 = rc.envClear;
                    if (Math.abs(pe0.getWidth() - cellW) > cellW*0.03 ||
                        Math.abs(pe0.getHeight()- cellH) > cellH*0.03) continue;

                    double dxa = anchorX - pe0.getMinX();
                    double dya = anchorY - pe0.getMinY();

                    Geometry ccTry = GeomUtils.translate(rc.clear, dxa, dya);
                    if (!sheetInset.covers(ccTry)) continue;

                    boolean clash = false;
                    for (Geometry pg : placedClear) if (collideStrict(pg, ccTry)) { clash = true; break; }
                    if (clash) continue;

                    Geometry cand = GeomUtils.translate(rc.base, dxa, dya);

                    PlacedPart pp = new PlacedPart();
                    pp.id = pi.spec.id; pp.polyPlaced = cand; pp.angleDeg = angDeg;
                    pp.mirrored = false; pp.dx = dxa; pp.dy = dya;

                    placed.add(pp);
                    placedClear.add(ccTry);
                    remaining.put(pi.spec.id, rem-1);
                    break;
                }
            }
        }
    }

    // =================== Utilidades varias ===================
    static double cv(double[] a){
        if (a==null || a.length==0) return 1.0;
        double s=0; for(double x:a) s+=x;
        double m=s/a.length; if (m==0) return 1.0;
        double v=0; for(double x:a) v+=(x-m)*(x-m);
        v/=a.length; return Math.sqrt(v)/m;
    }

    static double median(double[] a){
        if (a==null || a.length==0) return 0.0;
        double[] cp=a.clone(); Arrays.sort(cp); return cp[cp.length/2];
    }

    static void initLogger() {
        try {
            tee = new PrintStream(new FileOutputStream(LOG_FILE, false), true, StandardCharsets.UTF_8);
        } catch (Exception e) {
            tee = null;
        }
        log("[" + TS.format(new Date()) + "] Logger abierto");
    }
    static void closeLogger() { if (tee!=null) tee.close(); }
    static void log(String s) {
        String msg = TS.format(new Date()) + " | " + s;
        System.out.println(msg);
        if (tee!=null) tee.println(msg);
    }

    static Individual cloneShallow(Individual a){
        Individual b = new Individual();
        b.order = a.order.clone();
        b.angles = a.angles.clone();
        b.layout = a.layout;
        b.area = a.area; b.fitness = a.fitness; b.yMax = a.yMax;
        b.tieNoise = a.tieNoise;
        b.fingerprint = a.fingerprint;
        return b;
    }

    static Individual cloneDeep(Individual a){
        Individual b = new Individual();
        b.order = a.order.clone();
        b.angles = a.angles.clone();
        if (a.layout != null) b.layout = new ArrayList<>(a.layout);
        b.area = a.area; b.fitness = a.fitness; b.yMax = a.yMax;
        b.tieNoise = a.tieNoise;
        b.fingerprint = a.fingerprint;
        return b;
    }

    // =================== Fingerprint ===================
    static void registerFingerprint(Individual ind){
        long h = 1469598103934665603L; // FNV-like
        for (int v : ind.order)  { h ^= v; h *= 1099511628211L; }
        for (int v : ind.angles) { h ^= v*1315423911L; h *= 1099511628211L; }
        ind.fingerprint = h;
    }

    // =================== OX helpers ===================
    static int[] ox(int[] a, int[] b, Random rng){
        int n=a.length;
        int i = rng.nextInt(n), j = rng.nextInt(n);
        if (i>j){ int t=i; i=j; j=t; }
        int[] child = new int[n];
        Arrays.fill(child, -1);
        for (int k=i;k<=j;k++) child[k]=a[k];
        int cur=(j+1)%n;
        for (int k=0;k<n;k++) {
            int gene = b[(j+1+k)%n];
            if (!contains(child, gene)) { child[cur]=gene; cur=(cur+1)%n; }
        }
        return child;
    }

    static boolean[] originMaskOX(int[] p1, int[] p2, int[] child){
        boolean[] mask = new boolean[child.length];
        for (int i=0;i<child.length;i++) mask[i] = (child[i]==p1[i]);
        return mask;
    }

    static int indexOf(int[] arr, int val){
        for (int i=0;i<arr.length;i++) if (arr[i]==val) return i;
        return -1;
    }
    static boolean contains(int[] arr, int val){
        for (int v: arr) if (v==val) return true;
        return false;
    }
    static void swap(int[] a, int i, int j){ int t=a[i]; a[i]=a[j]; a[j]=t; }
    static void reverse(int[] a, int i, int j){ while (i<j){ swap(a,i++,j--); } }
    static void shuffle(int[] a, Random rng){
        for (int i=a.length-1;i>0;i--){
            int j=rng.nextInt(i+1);
            swap(a,i,j);
        }
    }
    static int norm(int ang){ int x = ang % 360; if (x<0) x+=360; return x; }

    static int snapToDomain(int ang, int[] domain){
        if (domain==null || domain.length==0) return ang;
        int best=domain[0], bestd = Math.abs(norm(ang)-domain[0]);
        for (int v: domain){
            int d = Math.abs(norm(ang)-v);
            if (d<bestd){ best=v; bestd=d; }
        }
        return best;
    }

    static List<PartInstance> expandInstances(List<PartSpec> specs) {
        List<PartInstance> out = new ArrayList<>();
        for (PartSpec s : specs) {
            int q = Math.max(1, s.qty);
            for (int i=0;i<q;i++) out.add(new PartInstance(s, i));
        }
        return out;
    }

    static Map<PartSpec,int[]> computeAngleOptions(List<PartSpec> specs) {
        Map<PartSpec,int[]> m = new HashMap<>();
        for (PartSpec s : specs) {
            int step = Math.max(1, s.rotationStepDeg);
            ArrayList<Integer> a = new ArrayList<>();
            for (int ang=0; ang<360; ang+=step) a.add(ang);
            int[] arr = new int[a.size()];
            for (int i=0;i<a.size();i++) arr[i] = a.get(i);
            m.put(s, arr);
        }
        return m;
    }

    static long randomSeed() {
        long s = System.nanoTime() ^ Double.doubleToLongBits(Math.random());
        return (s == 0L ? 1L : s);
    }
}
