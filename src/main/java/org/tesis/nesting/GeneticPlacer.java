package org.tesis.nesting;

import org.locationtech.jts.geom.Envelope;
import org.locationtech.jts.geom.Geometry;
import org.locationtech.jts.geom.Polygon;

import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class GeneticPlacer {

    // parámetros del algoritmo genético
    static final int    POP_SIZE            = 60;
    static final int    GENERATIONS_MAX     = 200;
    static final int    TOURNAMENT_K        = 3;
    static final double MUT_SWAP_P          = 0.35;
    static final double MUT_ANGLE_P         = 0.25;
    static final int    ELITE_COUNT         = 2;
    static final int    IMMIGRANTS_PER_GEN  = 2;
    static final int    EARLY_STOP_GENS     = 25;
    static final long   DEFAULT_SEED        = 42L;
    static final long   MAX_MILLIS_TOTAL    = 45_000;

    // parámetros del decodificador y colocación
    static final int    GRID_ANCHORS_X      = 10;
    static final int    GRID_ANCHORS_Y      = 6;
    static final int    MAX_ANCHORS_PER_TRY = 240;
    static final int    MAX_CHECKS_PER_PART = 1200;
    static final double NUDGE_STEP          = 0.35;
    static final int    NUDGE_STEPS_MAX     = 6;
    static final double CLEAR_SIMPLIFY_TOL  = 0.0;

    // pesos para evaluar posiciones
    static final double W_BEST_Y            = 1.0;
    static final double W_BEST_X            = 0.15;
    static final double W_BBOX_INFLATE      = 0.001;

    // peso poblacional
    static final double LAMBDA_HEIGHT       = 0.02;

    // modo lattice (rejilla)
    static final double GRID_CV_MAX         = 0.02;
    static final double SNAP_FRAC           = 0.03;
    static final double W_GRID_ALIGN        = 0.30;
    static final double W_CONTACT_LEFT      = 0.10;
    static final double W_CONTACT_DOWN      = 0.12;

    // contacto en modo no-rejilla
    static final double CONTACT_EPS         = 0.25;
    static final double W_TOUCH_LEFT_NG     = 0.90;
    static final double W_TOUCH_DOWN_NG     = 0.75;

    // compactación global
    static final int    POST_ROUNDS         = 2;
    static final int    POST_STEPS_MAX      = 8000;
    static final double POST_MIN_STEP       = 0.5;
    static final double POST_STEP_FRAC      = 0.08;

    // relleno de huecos
    static final int    FILL_ROUNDS         = 2;
    static final int    FILL_MAX_ANCHORS    = 2000;

    // logging
    static final String LOG_FILE = "LogGenetic.log";
    static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");
    static final DecimalFormat DF = new DecimalFormat("#,##0.###");
    private static PrintStream tee;

    // clases auxiliares
    static class PartInstance {
        final PartSpec spec; final int copyIdx;
        PartInstance(PartSpec s, int k) { spec = s; copyIdx = k; }
    }
    static class Individual {
        int[] order;
        int[] angleIdx;
        double fitnessArea;
        double usedYMax;
        double score;
        List<PlacedPart> placed;
    }
    static class Anchor {
        final double x, y;
        Anchor(double x, double y){ this.x=x; this.y=y; }
    }

    // API principal para colocar piezas con seed por defecto
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        return place(sheetInset, specs, DEFAULT_SEED);
    }

    // API principal para colocar piezas con seed definido
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs, long seed) {
        initLogger();
        final long t0 = System.currentTimeMillis();
        log("==== GeneticPlacer RUN " + new Date() + " ====");
        log("INI | PartSpec: " + specs.size() + " | seed=" + seed);

        // expandir copias de piezas
        List<PartInstance> instances = expandInstances(specs);
        if (instances.isEmpty()) { 
            log("No hay copias a colocar. Fin."); 
            closeLogger();
            return Collections.emptyList(); 
        }
        log("Copias expandidas: " + instances.size());

        // calcular ángulos disponibles por pieza
        Map<PartSpec,int[]> angleOptions = computeAngleOptions(specs);

        // crear población inicial
        Random rng = new Random(seed);
        List<Individual> pop = initPopulation(instances, angleOptions, rng);
        log("Población inicial: " + pop.size());

        // evaluación inicial
        evalPopulation(pop, sheetInset, instances, angleOptions);
        Individual best = bestOf(pop);
        log("Pre-loop | area=" + DF.format(best.fitnessArea) +
            " | yMax=" + DF.format(best.usedYMax) +
            " | score=" + DF.format(best.score) +
            " | colocadas=" + (best.placed==null?0:best.placed.size()));

        long lastImproveGen = -1;
        for (int gen=0; gen<GENERATIONS_MAX; gen++) {
            if (System.currentTimeMillis()-t0 > MAX_MILLIS_TOTAL) { log("Tiempo agotado. Gen " + gen); break; }
            long tg = System.currentTimeMillis();
            log("=== Gen " + gen + " INI ===");

            // elitismo (copiar los mejores)
            List<Individual> next = new ArrayList<>(POP_SIZE);
            for (Individual e : topK(pop, ELITE_COUNT)) next.add(cloneInd(e));

            // reproducción hasta llenar la población
            while (next.size() < POP_SIZE - IMMIGRANTS_PER_GEN) {
                Individual p1 = tournament(pop, rng);
                Individual p2 = tournament(pop, rng);
                Individual c  = crossover(p1, p2, rng);
                mutate(c, instances, angleOptions, rng);
                next.add(c);
            }
            // agregar inmigrantes aleatorios
            for (int i=0;i<IMMIGRANTS_PER_GEN;i++) next.add(randomIndividual(instances, angleOptions, rng));

            // evaluar nueva generación
            evalPopulation(next, sheetInset, instances, angleOptions);

            Individual genBest = bestOf(next);
            if (isBetter(genBest, best)) {
                best = cloneInd(genBest);
                lastImproveGen = gen;
                log("Mejora | area=" + DF.format(best.fitnessArea) +
                    " | yMax=" + DF.format(best.usedYMax) +
                    " | score=" + DF.format(best.score) +
                    " | colocadas=" + (best.placed==null?0:best.placed.size()));
            } else {
                log("Sin mejora | best score=" + DF.format(best.score));
            }

            if (lastImproveGen >= 0 && gen - lastImproveGen >= EARLY_STOP_GENS) {
                log("Early-stop tras " + EARLY_STOP_GENS + " generaciones sin mejora."); break;
            }
            pop = next;
            log("=== Gen " + gen + " FIN | dur=" + (System.currentTimeMillis()-tg) + " ms ===");
        }

        // resumen final (mismo formato que GRASP)
        int totalPiezas        = instances.size();
        int colocadas          = (best.placed == null ? 0 : best.placed.size());
        double areaPlancha     = sheetInset.getArea();
        double areaColocada    = Math.max(0.0, best.fitnessArea);
        double areaDesperdiciada = Math.max(0.0, areaPlancha - areaColocada);
        double aprovechamiento = areaPlancha > 0 ? (100.0 * areaColocada / areaPlancha) : 0.0;
        long durTotalMs        = System.currentTimeMillis() - t0;

        log("------------------------------");
        log("RESUMEN FINAL (GENETIC)");
        log("Seed               : " + seed);
        log("Piezas a posicionar: " + totalPiezas);
        log("Piezas colocadas   : " + colocadas);
        log("% aprovechamiento  : " + DF.format(aprovechamiento) + " %");
        log("Área colocada      : " + DF.format(areaColocada));
        log("Área plancha       : " + DF.format(areaPlancha));
        log("Área desperdiciada : " + DF.format(areaDesperdiciada));
        log("Ymax               : " + DF.format(best.usedYMax));
        log("Tiempo total       : " + durTotalMs + " ms (" + DF.format(durTotalMs/1000.0) + " s)");
        log("------------------------------");

        log("BEST | colocadas=" + colocadas +
            " | area=" + DF.format(areaColocada) +
            " | yMax=" + DF.format(best.usedYMax) +
            " | score=" + DF.format(best.score));

        closeLogger();
        return best.placed != null ? best.placed : Collections.emptyList();
    }

    // chequeo de colisión estricto (intersección real, no solo toque)
    static boolean collideStrict(Geometry a, Geometry b) {
        return a.intersects(b) && !a.touches(b);
    }

    // inicializa la población de individuos
    static List<Individual> initPopulation(List<PartInstance> instances,
                                           Map<PartSpec,int[]> angleOpts,
                                           Random rng) {
        List<Individual> pop = new ArrayList<>(POP_SIZE);
        for (int i=0;i<POP_SIZE;i++) pop.add(randomIndividual(instances, angleOpts, rng));

        // primer individuo greedy: orden por área descendente
        int n = instances.size();
        Integer[] idx = new Integer[n];
        for (int i=0;i<n;i++) idx[i]=i;
        Arrays.sort(idx, (a,b)->Double.compare(
                instances.get(b).spec.poly.getArea(),
                instances.get(a).spec.poly.getArea()
        ));
        Individual greedy = new Individual();
        greedy.order = new int[n];
        for (int i=0;i<n;i++) greedy.order[i] = idx[i];
        greedy.angleIdx = new int[n];
        for (int i=0;i<n;i++) {
            int[] ao = angleOpts.get(instances.get(i).spec);
            greedy.angleIdx[i] = (ao==null||ao.length==0)?0:rng.nextInt(ao.length);
        }
        pop.set(0, greedy);
        return pop;
    }

    // genera un individuo aleatorio
    static Individual randomIndividual(List<PartInstance> instances,
                                       Map<PartSpec,int[]> angleOpts,
                                       Random rng) {
        int n = instances.size();
        Individual ind = new Individual();
        ind.order = new int[n];
        for (int i=0;i<n;i++) ind.order[i]=i;
        for (int i=n-1;i>0;i--) { int j=rng.nextInt(i+1); int t=ind.order[i]; ind.order[i]=ind.order[j]; ind.order[j]=t; }
        ind.angleIdx = new int[n];
        for (int i=0;i<n;i++) {
            int[] ao = angleOpts.get(instances.get(i).spec);
            ind.angleIdx[i] = (ao==null||ao.length==0)?0:rng.nextInt(ao.length);
        }
        return ind;
    }

    // torneo de selección
    static Individual tournament(List<Individual> pop, Random rng) {
        Individual best = null;
        for (int i=0;i<TOURNAMENT_K;i++) {
            Individual c = pop.get(rng.nextInt(pop.size()));
            if (best==null || isBetter(c,best)) best=c;
        }
        return best;
    }

    // compara si un individuo es mejor que otro
    static boolean isBetter(Individual a, Individual b) {
        if (b==null) return true;
        return a.score > b.score;
    }

    // cruce entre dos individuos (OX para orden y 1 punto para ángulos)
    static Individual crossover(Individual p1, Individual p2, Random rng) {
        final int n = p1.order.length;
        Individual c = new Individual();

        int c1=rng.nextInt(n), c2=rng.nextInt(n);
        if (c1>c2){int t=c1; c1=c2; c2=t;}
        int[] child = new int[n];
        Arrays.fill(child,-1);
        boolean[] used = new boolean[n];

        for (int i=c1;i<=c2;i++){ int g=p1.order[i]; child[i]=g; used[g]=true; }

        int write=(c2+1)%n;
        for (int k=0;k<n;k++){
            int gene = p2.order[(c2+1+k)%n];
            if(!used[gene]){
                while(child[write]!=-1) write=(write+1)%n;
                child[write]=gene; used[gene]=true; write=(write+1)%n;
            }
        }
        c.order = child;

        int cut=rng.nextInt(n);
        c.angleIdx = new int[n];
        for (int i=0;i<n;i++) c.angleIdx[i] = (i<cut)? p1.angleIdx[i] : p2.angleIdx[i];
        return c;
    }

    // mutación de un individuo
    static void mutate(Individual ind,
                       List<PartInstance> instances,
                       Map<PartSpec,int[]> angleOpts,
                       Random rng) {
        int n = ind.order.length;

        if (rng.nextDouble()<MUT_SWAP_P && n>=2){
            int i=rng.nextInt(n), j=rng.nextInt(n);
            int t=ind.order[i]; ind.order[i]=ind.order[j]; ind.order[j]=t;
            int ta=ind.angleIdx[i]; ind.angleIdx[i]=ind.angleIdx[j]; ind.angleIdx[j]=ta;
        }
        for (int inst=0; inst<n; inst++){
            if (rng.nextDouble()<MUT_ANGLE_P){
                int[] ao = angleOpts.get(instances.get(inst).spec);
                ind.angleIdx[inst] = (ao==null||ao.length==0)?0:rng.nextInt(ao.length);
            }
        }
    }

    // evalúa la población
    static void evalPopulation(List<Individual> pop,
                               Polygon sheetInset,
                               List<PartInstance> instances,
                               Map<PartSpec,int[]> angleOptions) {
        int i=0;
        for (Individual ind: pop){
            log("Eval " + (++i) + "/" + pop.size());
            decode(ind, sheetInset, instances, angleOptions);
            ind.score = ind.fitnessArea - LAMBDA_HEIGHT * ind.usedYMax;
        }
    }

    // decodifica un individuo y construye su layout
    static void decode(Individual ind,
                       Polygon sheetInset,
                       List<PartInstance> instances,
                       Map<PartSpec,int[]> angleOptions) {
        long t0 = System.currentTimeMillis();

        final Envelope env = sheetInset.getEnvelopeInternal();
        final double minX = env.getMinX(), maxX = env.getMaxX();
        final double minY = env.getMinY(), maxY = env.getMaxY();

        // array para marcar qué piezas se colocaron en este pase
        final boolean[] wasPlaced = new boolean[ind.order.length];

        // detectar si aplica modo lattice (rejilla regular)
        boolean GRID_MODE = false;
        double cellW = 0, cellH = 0;
        int limit = Math.min(50, instances.size());
        if (limit > 0) {
            double[] ws = new double[limit], hs = new double[limit];
            for (int i=0;i<limit;i++){
                PartSpec sp = instances.get(i).spec;
                int[] allowed = angleOptions.get(sp);
                int angDeg = (allowed!=null && allowed.length>0) ? allowed[0] : 0;
                Geometry g0  = GeomUtils.orient(sp.poly, angDeg, false);
                Geometry gc0 = (sp.clearance > 0) ? g0.buffer(sp.clearance, 8) : g0;
                Envelope e = gc0.getEnvelopeInternal();
                ws[i] = e.getWidth(); hs[i] = e.getHeight();
            }
            if (cv(ws) <= GRID_CV_MAX && cv(hs) <= GRID_CV_MAX) {
                Envelope se = sheetInset.getEnvelopeInternal();
                double rectArea = se.getWidth()*se.getHeight();
                GRID_MODE = Math.abs(sheetInset.getArea() - rectArea) < 1e-6;
                if (GRID_MODE) { cellW = median(ws); cellH = median(hs); }
            }
        }
        final double snapEpsLocal = GRID_MODE ? Math.max(1e-6, SNAP_FRAC * Math.min(cellW, cellH)) : 0.0;

        // fast-path si estamos en modo rejilla
        if (GRID_MODE) {
            final List<PlacedPart> placed = new ArrayList<>();
            final List<Geometry>   placedClear = new ArrayList<>();

            gridFillRaster(
                sheetInset, placed, placedClear, instances, angleOptions,
                minX, minY, maxX, maxY, cellW, cellH
            );

            if (!placed.isEmpty()) postCompact(sheetInset, placed, placedClear);

            double totalArea = 0.0;
            double yMaxFinal = minY;
            for (PlacedPart pp : placed) {
                totalArea += pp.polyPlaced.getArea();
                yMaxFinal = Math.max(yMaxFinal, pp.polyPlaced.getEnvelopeInternal().getMaxY());
            }
            ind.placed      = placed;
            ind.fitnessArea = totalArea;
            ind.usedYMax    = yMaxFinal;

            log("  decode() GRID fast-path | colocadas=" + placed.size() +
                " | área=" + DF.format(ind.fitnessArea) +
                " | yMax=" + DF.format(ind.usedYMax) +
                " | dur=" + (System.currentTimeMillis()-t0) + " ms");
            return;
        }

        // flujo normal (cuando no aplica rejilla)
        final List<PlacedPart> placed = new ArrayList<>();
        final List<Geometry>   placedClear = new ArrayList<>();
        final ArrayList<Anchor> anchors = initialAnchors(minX, minY, maxX, maxY);

        double globMinX = maxX, globMinY = maxY, globMaxX = minX, globMaxY = minY;
        double areaSum = 0.0;
        double usedYMax = minY;

        Random rndPick = new Random(ind.hashCode() ^ 0x51F15E5D);

        // recorre cada pieza en el orden del individuo
        for (int pos=0; pos<ind.order.length; pos++) {
            int instIdx = ind.order[pos];
            PartSpec sp = instances.get(instIdx).spec;

            int[] allowed = angleOptions.get(sp);
            int angDeg = (allowed!=null && allowed.length>0) ? allowed[ind.angleIdx[instIdx] % allowed.length] : 0;

            Geometry g0  = GeomUtils.orient(sp.poly, angDeg, false);
            Geometry g   = g0;
            Geometry gc0 = (sp.clearance > 0) ? g0.buffer(sp.clearance, 8) : g0;
            Geometry gc  = (CLEAR_SIMPLIFY_TOL>0) ? gc0.buffer(0, 0) : gc0;

            Envelope pe = gc.getEnvelopeInternal();
            double width = pe.getWidth(), height = pe.getHeight();
            if (width <= 0 || height <= 0) continue;

            double x0 = minX, x1 = maxX - width;
            double y0 = minY, y1 = maxY - height;
            if (x1 < x0 || y1 < y0) continue;

            // elegir subconjunto de anclas a probar
            List<Anchor> pool = sampleAnchors(anchors, MAX_ANCHORS_PER_TRY, rndPick);

            int checks = 0;
            double bestVal = Double.POSITIVE_INFINITY;
            Geometry bestCand=null, bestClear=null;
            double bestDX=0, bestDY=0;

            // probar cada ancla
            for (Anchor a : pool) {
                if (checks > MAX_CHECKS_PER_PART) break;

                double dxa = a.x - pe.getMinX();
                double dya = a.y - pe.getMinY();

                double txMin = pe.getMinX() + dxa, txMax = pe.getMaxX() + dxa;
                double tyMin = pe.getMinY() + dya, tyMax = pe.getMaxY() + dya;
                if (txMin < minX-1e-9 || txMax > maxX+1e-9 || tyMin < minY-1e-9 || tyMax > maxY+1e-9) continue;

                Geometry cc = GeomUtils.translate(gc, dxa, dya);
                checks++;
                if (!sheetInset.covers(cc)) continue;

                boolean inter = false;
                for (Geometry pg : placedClear) {
                    if (collideStrict(pg, cc)) { inter = true; break; }
                }
                if (inter) continue;

                // intentar empujar la pieza hacia la izquierda
                double[] nudged = tryNudge(sheetInset, cc, placedClear, -NUDGE_STEP, 0, NUDGE_STEPS_MAX);
                cc = GeomUtils.translate(cc, nudged[0], nudged[1]);
                Geometry cand = GeomUtils.translate(g,  dxa+nudged[0], dya+nudged[1]);

                // cálculo de score local
                Envelope ce = cc.getEnvelopeInternal();
                double cy = ce.getMinY();
                double cx = ce.getMinX();

                double nGlobMinX = Math.min(globMinX, ce.getMinX());
                double nGlobMinY = Math.min(globMinY, ce.getMinY());
                double nGlobMaxX = Math.max(globMaxX, ce.getMaxX());
                double nGlobMaxY = Math.max(globMaxY, ce.getMaxY());
                double bboxGrow = (nGlobMaxX - nGlobMinX)*(nGlobMaxY - nGlobMinY)
                                - Math.max(0, (globMaxX - globMinX)*(globMaxY - globMinY));

                // bonos por contacto
                double tol = CONTACT_EPS;
                double touchL = 0, touchD = 0;

                if (Math.abs(ce.getMinX() - minX) <= tol) touchL = 1;
                if (Math.abs(ce.getMinY() - minY) <= tol) touchD = 1;

                for (Geometry pg : placedClear) {
                    Envelope pe2 = pg.getEnvelopeInternal();
                    if (Math.abs(ce.getMinX() - pe2.getMaxX()) <= tol &&
                       !(ce.getMaxY() <= pe2.getMinY() || ce.getMinY() >= pe2.getMaxY()))
                        touchL = 1;
                    if (Math.abs(ce.getMinY() - pe2.getMaxY()) <= tol &&
                       !(ce.getMaxX() <= pe2.getMinX() || ce.getMinX() >= pe2.getMaxX()))
                        touchD = 1;
                    if (touchL==1 && touchD==1) break;
                }
                double contactBonus = - (W_TOUCH_LEFT_NG * touchL + W_TOUCH_DOWN_NG * touchD);

                double val = W_BEST_Y*cy + W_BEST_X*cx + W_BBOX_INFLATE*bboxGrow + contactBonus;

                if (val < bestVal) {
                    bestVal = val;
                    bestCand = cand;
                    bestClear= cc;
                    bestDX = dxa+nudged[0];
                    bestDY = dya+nudged[1];
                }
            }

            // si encontramos una posición válida, fijamos la pieza
            if (bestCand != null) {
                PlacedPart pp = new PlacedPart();
                pp.id = sp.id;
                pp.polyPlaced = bestCand;
                pp.angleDeg = angDeg;
                pp.mirrored = false;
                pp.dx = bestDX; pp.dy = bestDY;

                placed.add(pp);
                placedClear.add(bestClear);
                wasPlaced[pos] = true;

                areaSum += bestCand.getArea();
                usedYMax = Math.max(usedYMax, bestCand.getEnvelopeInternal().getMaxY());

                Envelope ce = bestClear.getEnvelopeInternal();
                globMinX = Math.min(globMinX, ce.getMinX());
                globMinY = Math.min(globMinY, ce.getMinY());
                globMaxX = Math.max(globMaxX, ce.getMaxX());
                globMaxY = Math.max(globMaxY, ce.getMaxY());

                addLocalAnchors(anchors, ce, minX, minY, maxX, maxY);
            }

            if ((pos % 10) == 0) {
                log("    · pieza " + pos + "/" + ind.order.length +
                        " | angle=" + angDeg +
                        " | areaAcum=" + DF.format(areaSum) +
                        " | anchors=" + anchors.size());
            }
        }

        // aplicar compactación global
        postCompact(sheetInset, placed, placedClear);

        // intentar rellenar huecos con piezas no colocadas
        ArrayList<PartInstance> unplaced = new ArrayList<>();
        for (int k = 0; k < ind.order.length; k++) {
            if (!wasPlaced[k]) unplaced.add(instances.get(ind.order[k]));
        }
        if (!unplaced.isEmpty()) {
            int added = fillGapsGreedy(sheetInset, placed, placedClear, unplaced, angleOptions);
            if (added > 0) postCompact(sheetInset, placed, placedClear);
        }

        // métricas finales del individuo
        double totalArea = 0.0;
        double yMaxFinal = minY;
        for (PlacedPart pp : placed) {
            totalArea += pp.polyPlaced.getArea();
            yMaxFinal = Math.max(yMaxFinal, pp.polyPlaced.getEnvelopeInternal().getMaxY());
        }

        ind.placed      = placed;
        ind.fitnessArea = totalArea;
        ind.usedYMax    = yMaxFinal;

        log("  decode() end | colocadas=" + placed.size() +
            " | área=" + DF.format(ind.fitnessArea) +
            " | yMax=" + DF.format(ind.usedYMax) +
            " | dur=" + (System.currentTimeMillis()-t0) + " ms");
    }

    // crea anclas iniciales (esquinas + malla)
    static ArrayList<Anchor> initialAnchors(double minX, double minY, double maxX, double maxY){
        ArrayList<Anchor> a = new ArrayList<>();
        a.add(new Anchor(minX, minY));
        a.add(new Anchor(maxX, minY));
        a.add(new Anchor(minX, maxY));
        a.add(new Anchor(maxX, maxY));

        for (int i=0;i<=GRID_ANCHORS_X;i++){
            double fx = i/(double)GRID_ANCHORS_X;
            double x = minX + fx*(maxX-minX);
            for (int j=0;j<=GRID_ANCHORS_Y;j++){
                double fy = j/(double)GRID_ANCHORS_Y;
                double y = minY + fy*(maxY-minY);
                a.add(new Anchor(x,y));
            }
        }
        return a;
    }

    // añade anclas locales alrededor de una pieza colocada
    static void addLocalAnchors(List<Anchor> anchors, Envelope ce, double minX, double minY, double maxX, double maxY){
        double left = ce.getMinX(), right = ce.getMaxX();
        double bottom = ce.getMinY(), top = ce.getMaxY();

        if (right <= maxX) {
            anchors.add(new Anchor(right, bottom));
            anchors.add(new Anchor(right, top));
            anchors.add(new Anchor(right, (bottom+top)*0.5));
        }
        if (top <= maxY) {
            anchors.add(new Anchor(left,  top));
            anchors.add(new Anchor(right, top));
            anchors.add(new Anchor((left+right)*0.5, top));
        }
        double eps = 1e-6;
        anchors.add(new Anchor(Math.min(right, maxX)-eps, bottom));
        anchors.add(new Anchor(left, Math.min(top, maxY)-eps));
    }

    // devuelve un subconjunto aleatorio de anclas
    static List<Anchor> sampleAnchors(List<Anchor> all, int limit, Random rng){
        if (all.size() <= limit) return new ArrayList<>(all);
        ArrayList<Anchor> copy = new ArrayList<>(all);
        for (int i=copy.size()-1; i>0; i--){
            int j = rng.nextInt(i+1);
            Anchor t = copy.get(i); copy.set(i, copy.get(j)); copy.set(j, t);
        }
        return copy.subList(0, limit);
    }

    // intenta desplazar una geometría paso a paso (nudge)
    static double[] tryNudge(Polygon sheet, Geometry clearGeom, List<Geometry> placedClear,
                             double stepX, double stepY, int stepsMax) {
        double accX = 0, accY = 0;
        Geometry g = clearGeom;
        for (int s=0; s<stepsMax; s++){
            double nx = accX + stepX;
            double ny = accY + stepY;
            Geometry t = GeomUtils.translate(clearGeom, nx, ny);
            if (!sheet.contains(t)) break;
            boolean inter=false;
            for (Geometry pg : placedClear){ if (pg.intersects(t)){ inter=true; break; } }
            if (inter) break;
            accX = nx; accY = ny; g = t;
        }
        if (stepX!=0 && stepY!=0){
            double[] down = tryNudge(sheet, g, placedClear, 0, stepY, stepsMax/2);
            accX += down[0]; accY += down[1];
        }
        return new double[]{accX, accY};
    }

    // calcula los ángulos posibles de cada pieza
    static Map<PartSpec,int[]> computeAngleOptions(List<PartSpec> specs) {
        Map<PartSpec,int[]> map = new HashMap<>();
        for (PartSpec sp : specs) {
            int step = Math.max(0, sp.rotationStepDeg);
            int[] arr;
            if (step <= 0) {
                arr = new int[]{0};
            } else {
                int cnt = 360 / gcd(360, step);
                arr = new int[cnt];
                int a=0; for (int i=0;i<cnt;i++){ arr[i]=a; a=(a+step)%360; }
            }
            map.put(sp, arr);
        }
        return map;
    }

    // expande la lista de PartSpec a instancias individuales
    static List<PartInstance> expandInstances(List<PartSpec> specs) {
        List<PartInstance> out = new ArrayList<>();
        for (PartSpec sp : specs) {
            int copies = Math.max(1, sp.qty);
            for (int k=0;k<copies;k++) out.add(new PartInstance(sp,k));
        }
        return out;
    }
    // devuelve el mejor individuo de la población
    static Individual bestOf(List<Individual> pop){
        Individual best=null; for (Individual ind:pop) if (isBetter(ind,best)) best=ind; return best;
    }

    // devuelve el top K de la población
    static List<Individual> topK(List<Individual> pop, int k){
        ArrayList<Individual> cp=new ArrayList<>(pop);
        cp.sort((a,b)->Double.compare(b.score,a.score));
        return cp.subList(0, Math.min(k, cp.size()));
    }

    // clona un individuo
    static Individual cloneInd(Individual a){
        Individual b=new Individual();
        b.order=a.order.clone();
        b.angleIdx=a.angleIdx.clone();
        b.fitnessArea=a.fitnessArea;
        b.usedYMax=a.usedYMax;
        b.score=a.score;
        b.placed=a.placed;
        return b;
    }

    // máximo común divisor (para ángulos)
    static int gcd(int a,int b){ while(b!=0){ int t=a%b; a=b; b=t; } return Math.abs(a); }

    // helpers nuevos: grid snap, coeficiente de variación, mediana, etc.
    static double snapToGrid(double v, double origin, double step) {
        double t = (v - origin) / step;
        return origin + Math.rint(t) * step;
    }

    static double cv(double[] a){
        if (a==null || a.length==0) return 1.0;
        double s=0; for(double x:a) s+=x;
        double m=s/a.length; if (m==0) return 1.0;
        double v=0; for(double x:a) v+=(x-m)*(x-m);
        v/=a.length;
        return Math.sqrt(v)/m;
    }

    static double median(double[] a){
        if (a==null || a.length==0) return 0.0;
        double[] cp = a.clone();
        Arrays.sort(cp);
        return cp[cp.length/2];
    }

    // genera anclas en rejilla (modo lattice)
    static ArrayList<Anchor> latticeAnchors(double minX, double minY, double maxX, double maxY,
                                            double cellW, double cellH) {
        ArrayList<Anchor> a = new ArrayList<>();
        if (cellW <= 0 || cellH <= 0) return a;
        int nx = (int)Math.floor((maxX - minX) / cellW);
        int ny = (int)Math.floor((maxY - minY) / cellH);
        for (int j = 0; j < ny; j++) {
            double y = minY + j * cellH;
            for (int i = 0; i < nx; i++) {
                double x = minX + i * cellW;
                a.add(new Anchor(x, y));
            }
        }
        return a;
    }

    static long cellKey(int col, int row){ return (((long)col)<<32) ^ (row & 0xffffffffL); }
    static int colFromX(double x, double minX, double cellW){
        return (int)Math.floor((x - minX) / cellW + 1e-9);
    }
    static int rowFromY(double y, double minY, double cellH){
        return (int)Math.floor((y - minY) / cellH + 1e-9);
    }

    // chequeo estricto de solapamiento de envelopes
    static boolean envelopeOverlapsStrict(Envelope a, Envelope b, double eps) {
        return (a.getMinX() < b.getMaxX() - eps) &&
               (a.getMaxX() > b.getMinX() + eps) &&
               (a.getMinY() < b.getMaxY() - eps) &&
               (a.getMaxY() > b.getMinY() + eps);
    }

    // compactación global: intenta deslizar piezas a izquierda/abajo
    static double[] slide1D(Polygon sheet, Geometry g, List<Geometry> all, int idx,
                            double dirX, double dirY) {
        Envelope e = g.getEnvelopeInternal();
        double step = Math.max(POST_MIN_STEP, POST_STEP_FRAC * Math.min(e.getWidth(), e.getHeight()));
        double accX = 0, accY = 0;
        for (int s=0; s<POST_STEPS_MAX; s++){
            double nx = accX + dirX * step;
            double ny = accY + dirY * step;
            Geometry t = GeomUtils.translate(g, nx, ny);
            if (!sheet.contains(t)) break;
            boolean inter=false;
            for (int j=0; j<all.size(); j++){
                if (j==idx) continue;
                if (all.get(j).intersects(t)) { inter=true; break; }
            }
            if (inter) break;
            accX = nx; accY = ny;
        }
        return new double[]{accX, accY};
    }

    static void postCompact(Polygon sheet, List<PlacedPart> placed, List<Geometry> placedClear){
        for (int round=0; round<POST_ROUNDS; round++){
            // izquierda
            for (int i=0; i<placed.size(); i++){
                Geometry gi = placedClear.get(i);
                double[] d = slide1D(sheet, gi, placedClear, i, -1, 0);
                if (d[0]!=0 || d[1]!=0){
                    placedClear.set(i, GeomUtils.translate(gi, d[0], d[1]));
                    placed.get(i).polyPlaced = GeomUtils.translate(placed.get(i).polyPlaced, d[0], d[1]);
                    placed.get(i).dx += d[0]; placed.get(i).dy += d[1];
                }
            }
            // abajo
            for (int i=0; i<placed.size(); i++){
                Geometry gi = placedClear.get(i);
                double[] d = slide1D(sheet, gi, placedClear, i, 0, -1);
                if (d[0]!=0 || d[1]!=0){
                    placedClear.set(i, GeomUtils.translate(gi, d[0], d[1]));
                    placed.get(i).polyPlaced = GeomUtils.translate(placed.get(i).polyPlaced, d[0], d[1]);
                    placed.get(i).dx += d[0]; placed.get(i).dy += d[1];
                }
            }
        }
    }

    // intenta rellenar huecos con piezas que quedaron sin colocar
    static int fillGapsGreedy(Polygon sheetInset,
                              List<PlacedPart> placed,
                              List<Geometry> placedClear,
                              List<PartInstance> unplaced,
                              Map<PartSpec,int[]> angleOptions)
    {
        if (unplaced.isEmpty()) return 0;

        unplaced.sort((a,b) -> Double.compare(b.spec.poly.getArea(), a.spec.poly.getArea()));

        double[] ws = new double[unplaced.size()];
        double[] hs = new double[unplaced.size()];
        for (int i=0;i<unplaced.size();i++){
            Envelope e = unplaced.get(i).spec.poly.getEnvelopeInternal();
            ws[i] = e.getWidth(); hs[i] = e.getHeight();
        }
        double medW = Math.max(POST_MIN_STEP, 0.5 * Math.max(1e-9, median(ws)));
        double medH = Math.max(POST_MIN_STEP, 0.5 * Math.max(1e-9, median(hs)));

        Envelope env = sheetInset.getEnvelopeInternal();
        double minX = env.getMinX(), maxX = env.getMaxX();
        double minY = env.getMinY(), maxY = env.getMaxY();

        int placedNow = 0;

        for (int round = 0; round < FILL_ROUNDS && !unplaced.isEmpty(); round++) {
            ArrayList<Anchor> anchors = new ArrayList<>();

            int nx = (int)Math.max(1, Math.floor((maxX - minX) / medW));
            int ny = (int)Math.max(1, Math.floor((maxY - minY) / medH));
            for (int i=0; i<=nx; i++){
                double x = minX + i * medW;
                for (int j=0; j<=ny; j++){
                    double y = minY + j * medH;
                    anchors.add(new Anchor(x,y));
                }
            }

            for (Geometry g : placedClear) {
                Envelope ce = g.getEnvelopeInternal();
                double l = ce.getMinX(), r = ce.getMaxX();
                double b = ce.getMinY(), t = ce.getMaxY();
                anchors.add(new Anchor(r, b)); anchors.add(new Anchor(r, t));
                anchors.add(new Anchor((l+r)*0.5, t)); anchors.add(new Anchor(r, (b+t)*0.5));
                double eps = 1e-6;
                anchors.add(new Anchor(Math.min(r, maxX)-eps, b));
                anchors.add(new Anchor(l, Math.min(t, maxY)-eps));
            }

            if (anchors.size() > FILL_MAX_ANCHORS) {
                Collections.shuffle(anchors, new Random(1337+round));
                anchors.subList(FILL_MAX_ANCHORS, anchors.size()).clear();
            }

            Iterator<PartInstance> it = unplaced.iterator();
            while (it.hasNext()){
                PartInstance pi = it.next();
                PartSpec sp = pi.spec;

                int[] allowed = angleOptions.getOrDefault(sp, new int[]{0});
                boolean placedOk = false;

                for (int angDeg : allowed) {
                    Geometry g0  = GeomUtils.orient(sp.poly, angDeg, false);
                    Geometry gc0 = (sp.clearance > 0) ? g0.buffer(sp.clearance, 8) : g0;
                    Envelope pe0 = gc0.getEnvelopeInternal();

                    for (Anchor a : anchors) {
                        double dxa = a.x - pe0.getMinX();
                        double dya = a.y - pe0.getMinY();

                        Geometry ccTry = GeomUtils.translate(gc0, dxa, dya);
                        Envelope ce = ccTry.getEnvelopeInternal();

                        if (ce.getMinX() < minX-1e-9 || ce.getMaxX() > maxX+1e-9 ||
                            ce.getMinY() < minY-1e-9 || ce.getMaxY() > maxY+1e-9) continue;
                        if (!sheetInset.covers(ccTry)) continue;
                        boolean clash=false; for (Geometry pg: placedClear){ if (pg.intersects(ccTry)){ clash=true; break; } }
                        if (clash) continue;

                        double[] nudged = tryNudge(sheetInset, ccTry, placedClear, 0, -NUDGE_STEP, NUDGE_STEPS_MAX);
                        ccTry = GeomUtils.translate(ccTry, nudged[0], nudged[1]);

                        Geometry cand = GeomUtils.translate(g0, dxa+nudged[0], dya+nudged[1]);
                        PlacedPart pp = new PlacedPart();
                        pp.id = sp.id; pp.polyPlaced = cand; pp.angleDeg = angDeg; pp.mirrored=false; pp.dx=dxa+nudged[0]; pp.dy=dya+nudged[1];

                        placed.add(pp);
                        placedClear.add(ccTry);
                        placedNow++;
                        placedOk = true;
                        break;
                    }
                    if (placedOk) break;
                }

                if (placedOk) it.remove();
            }
        }

        return placedNow;
    }

    // rellenado en modo rejilla
    static void gridFillRaster(Polygon sheetInset,
                           List<PlacedPart> placed,
                           List<Geometry> placedClear,
                           List<PartInstance> instances,
                           Map<PartSpec,int[]> angleOptions,
                           double minX,double minY,double maxX,double maxY,
                           double cellW,double cellH)
    {
        final double EPS = 1e-9;
        int nx = Math.max(1, (int)Math.floor(((maxX - minX) + EPS) / cellW));
        int ny = Math.max(1, (int)Math.floor(((maxY - minY) + EPS) / cellH));

        Map<String,Integer> remaining = new HashMap<>();
        for (PartInstance pi : instances) remaining.merge(pi.spec.id, 1, Integer::sum);
        for (PlacedPart p : placed)       remaining.merge(p.id, -1, Integer::sum);

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
                        if (!sheetInset.covers(ccTry)) continue;

                        boolean clash = false;
                        for (Geometry pg : placedClear) {
                            if (collideStrict(pg, ccTry)) { clash = true; break; }
                        }
                        if (clash) continue;

                        Geometry cand = GeomUtils.translate(g0, dxa, dya);

                        PlacedPart pp = new PlacedPart();
                        pp.id = pi.spec.id; pp.polyPlaced = cand; pp.angleDeg = angDeg;
                        pp.mirrored = false; pp.dx = dxa; pp.dy = dya;

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

    // logger para guardar salida en archivo y consola
    static void initLogger(){
        if (tee!=null) return;
        try {
            PrintStream f = new PrintStream(new FileOutputStream(LOG_FILE, true), true, StandardCharsets.UTF_8);
            tee = new TeePrintStream(System.err, f);
            System.setOut(tee);
            System.setErr(tee);
            log("==== Logger inicializado (" + LOG_FILE + ") ====");
        } catch (Exception e) { e.printStackTrace(); }
    }

    static void closeLogger() {
        if (tee != null) tee.close();
    }

    static void log(String s){
        String line="["+TS.format(new Date())+"] "+s;
        (tee!=null?tee:System.err).println(line);
    }

    // clase para redirigir logs a dos streams
    static class TeePrintStream extends PrintStream {
        private final PrintStream a,b;
        TeePrintStream(PrintStream a, PrintStream b){ super(a); this.a=a; this.b=b; }
        @Override public void println(String x){ a.println(x); b.println(x); }
        @Override public void print(String s){ a.print(s); b.print(s); }
        @Override public void flush(){ a.flush(); b.flush(); }
        @Override public void close(){ try{a.close();}finally{b.close();} }
    }
}
