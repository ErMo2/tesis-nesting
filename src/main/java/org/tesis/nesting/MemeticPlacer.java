package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public class MemeticPlacer {

    static final int    POP_SIZE          = 36;
    static final int    ELITE_SIZE        = 4;
    static final int    GEN_MAX           = 200;
    static final long   MAX_MILLIS_TOTAL  = 60_000;
    static final long   MAX_MILLIS_LS     = 20_000;
    static final double MUT_RATE_ORDER    = 0.35;
    static final double MUT_RATE_ANGLE    = 0.25;
    static final int    TOURNAMENT_K      = 3;

    static final long   DEFAULT_SEED      = 42L;

    static final int GRID_ANCHORS_X      = 10;
    static final int GRID_ANCHORS_Y      = 6;
    static final int MAX_ANCHORS_PER_TRY = 220;
    static final int MAX_CHECKS_PER_PART = 1200;
    static final double NUDGE_STEP       = 0.5;
    static final int NUDGE_STEPS_MAX     = 6;

    static final double W_BEST_Y       = 1.0;
    static final double W_BEST_X       = 0.15;
    static final double W_BBOX_INFLATE = 0.01;

    static final double GRID_CV_MAX = 0.02;
    static final double EPS_EDGE    = 1e-9;

    static final String LOG_FILE = "LogMemetic.log";
    static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");
    static final DecimalFormat DF = new DecimalFormat("#,##0.###");
    private static PrintStream tee;

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
        double area;
        double yMax;
        long   tieNoise;
    }

    // ejecutar el algoritmo con semilla aleatoria para obtener resultados distintos en cada corrida.
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        long seed = randomSeed();
        return place(sheetInset, specs, seed);
    }

    // ejecutar el algoritmo con una semilla fija para reproducibilidad de resultados.
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs, long seed) {
        initLogger();
        final long t0 = System.currentTimeMillis();

        List<PartInstance> instances = expandInstances(specs);
        final int N = instances.size();
        if (N == 0) { log("No hay instancias"); closeLogger(); return Collections.emptyList(); }

        Map<PartSpec,int[]> angleOptionsBySpec = computeAngleOptions(specs);

        int[][] angleDomain = new int[N][];
        for (int i=0;i<N;i++) {
            angleDomain[i] = angleOptionsBySpec.get(instances.get(i).spec);
            if (angleDomain[i] == null || angleDomain[i].length == 0) angleDomain[i] = new int[]{0};
        }

        Random rng = new Random(seed);

        List<Individual> pop = new ArrayList<>(POP_SIZE);

        int greedyVariants = Math.min(ELITE_SIZE, POP_SIZE/3);
        for (int g=0; g<greedyVariants; g++) {
            pop.add(greedySeed(N, instances, angleDomain, rng));
        }
        while (pop.size() < POP_SIZE) {
            Individual ind = randomIndividual(N, angleDomain, rng);
            if (rng.nextDouble() < 0.5) {
                int i = rng.nextInt(ind.order.length-1);
                reverse(ind.order, i, Math.min(ind.order.length-1, i + 1 + rng.nextInt(Math.max(1, ind.order.length/5))));
            }
            pop.add(ind);
        }

        for (Individual ind : pop) evaluate(ind, sheetInset, instances, angleDomain, rng);
        Individual best = bestOf(pop);

        log("INI | seed=" + seed + " | N=" + N + " | pop=" + POP_SIZE + " | area=" + DF.format(best.area) + " | yMax=" + DF.format(best.yMax));

        for (int gen = 0; gen < GEN_MAX; gen++) {
            if (System.currentTimeMillis() - t0 > MAX_MILLIS_TOTAL) {
                log("Tiempo total agotado en gen=" + gen); break;
            }

            pop.sort(Comparator.<Individual>comparingDouble((Individual a)->-a.area)
                               .thenComparingDouble(a -> a.yMax)
                               .thenComparingLong(a -> a.tieNoise));
            List<Individual> next = new ArrayList<>(POP_SIZE);
            for (int e=0; e<ELITE_SIZE; e++) next.add(cloneShallow(pop.get(e)));

            while (next.size() < POP_SIZE) {
                Individual p1 = tournament(pop, rng, TOURNAMENT_K);
                Individual p2 = tournament(pop, rng, TOURNAMENT_K);
                Individual child = crossoverOX_withAngleBlend(p1, p2, angleDomain, rng);

                double mutOrder = MUT_RATE_ORDER * (0.85 + 0.3*rng.nextDouble());
                double mutAngle = MUT_RATE_ANGLE * (0.85 + 0.3*rng.nextDouble());
                mutate(child, angleDomain, rng, mutOrder, mutAngle);

                evaluate(child, sheetInset, instances, angleDomain, rng);
                next.add(child);
            }

            pop = next;

            long lsBudget = Math.max(0, MAX_MILLIS_LS - (System.currentTimeMillis() - t0) / 6);
            if (lsBudget > 0) {
                for (int e=0; e<Math.min(ELITE_SIZE, pop.size()); e++) {
                    localSearch(pop.get(e), sheetInset, instances, angleDomain, rng, lsBudget / ELITE_SIZE);
                }
            }

            if (gen>0 && gen % 25 == 0) {
                pop.sort(Comparator.<Individual>comparingDouble((Individual a)->a.area)
                                   .thenComparingDouble(a -> -a.yMax));
                int pert = Math.max(2, POP_SIZE/6);
                for (int i=0; i<pert; i++) {
                    Individual w = pop.get(i);
                    int len = Math.max(2, w.order.length/8);
                    int s = rng.nextInt(Math.max(1, w.order.length-len));
                    reverse(w.order, s, s+len-1);
                    for (int k=0; k<w.angles.length; k++) {
                        if (rng.nextDouble() < 0.15) {
                            int[] dom = angleDomain[w.order[k]];
                            w.angles[k] = dom[rng.nextInt(dom.length)];
                        }
                    }
                    evaluate(w, sheetInset, instances, angleDomain, rng);
                }
            }

            Individual curBest = bestOf(pop);
            if (isBetter(curBest, best)) best = cloneDeep(curBest);

            log(String.format("GEN %d | best area=%s yMax=%s | placed=%d",
                    gen, DF.format(best.area), DF.format(best.yMax),
                    (best.layout==null?0:best.layout.size())));
        }

        int totalPiezas        = instances.size();
        int colocadas          = (best.layout == null ? 0 : best.layout.size());
        double areaPlancha     = sheetInset.getArea();
        double areaColocada    = best.area;
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

    //decodificar una permutación y ángulos en un layout, usando fast-path en rejilla cuando aplica.
    static DecodeResult decode(Polygon sheetInset,
                               List<PartInstance> orderInst,
                               int[] angles,
                               Random rng)
    {
        DecodeResult res = new DecodeResult();
        Envelope env = sheetInset.getEnvelopeInternal();
        double minX = env.getMinX(), minY = env.getMinY(), maxX = env.getMaxX(), maxY = env.getMaxY();
        res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;

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
                    Geometry g0  = GeomUtils.orient(sp.poly, angDeg, false);
                    Geometry gc0 = (sp.clearance > 0) ? g0.buffer(sp.clearance, 8) : g0;
                    Envelope e   = gc0.getEnvelopeInternal();
                    ws[i] = e.getWidth(); hs[i] = e.getHeight();
                }
                if (cv(ws) <= GRID_CV_MAX && cv(hs) <= GRID_CV_MAX) {
                    GRID_MODE = true; cellW = median(ws); cellH = median(hs);
                }
            }
        }

        if (GRID_MODE) {
            gridFillRaster(sheetInset, res.placed, res.placedClear, orderInst, angles,
                           minX, minY, maxX, maxY, cellW, cellH);

            res.globMinX = maxX; res.globMinY = maxY; res.globMaxX = minX; res.globMaxY = minY;
            for (Geometry cc : res.placedClear) {
                Envelope ce = cc.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());
            }
            return res;
        }

        ArrayList<Anchor> anchors = initialAnchors(minX, minY, maxX, maxY);

        for (int pos = 0; pos < orderInst.size(); pos++) {
            PartInstance inst = orderInst.get(pos);
            int angDeg = norm(angles[pos]);

            Geometry gBase = GeomUtils.orient(inst.spec.poly, angDeg, false);
            Geometry gClr  = inst.spec.clearance > 0 ? gBase.buffer(inst.spec.clearance, 8) : gBase;
            Envelope pe    = gClr.getEnvelopeInternal();

            List<Anchor> sampleAnchors = sampleAnchors(anchors, rng, MAX_ANCHORS_PER_TRY);
            Candidate best = null;
            int checks=0;

            for (Anchor a : sampleAnchors) {
                if (++checks > MAX_CHECKS_PER_PART) break;

                double dxa = a.x - pe.getMinX();
                double dya = a.y - pe.getMinY();

                double txMin = pe.getMinX() + dxa, txMax = pe.getMaxX() + dxa;
                double tyMin = pe.getMinY() + dya, tyMax = pe.getMaxY() + dya;
                if (txMin < minX-EPS_EDGE || txMax > maxX+EPS_EDGE || tyMin < minY-EPS_EDGE || tyMax > maxY+EPS_EDGE) continue;

                Geometry cc = GeomUtils.translate(gClr, dxa, dya);
                if (!sheetInset.covers(cc)) continue;

                boolean inter = false;
                for (Geometry pg : res.placedClear) { if (collideStrict(pg, cc)) { inter = true; break; } }
                if (inter) continue;

                double[] nudged = tryNudge(sheetInset, cc, res.placedClear, -NUDGE_STEP, -NUDGE_STEP, NUDGE_STEPS_MAX);
                cc = GeomUtils.translate(cc, nudged[0], nudged[1]);
                Geometry placed = GeomUtils.translate(gBase, dxa+nudged[0], dya+nudged[1]);

                double val = localScore(res, cc);

                if (best == null || val < best.scoreVal) {
                    best = new Candidate();
                    best.clearGeom=cc; best.placedGeom=placed;
                    best.dx = dxa+nudged[0]; best.dy = dya+nudged[1];
                    best.angleDeg = angDeg; best.scoreVal = val;
                }
            }

            if (best != null) {
                PlacedPart pp = new PlacedPart();
                pp.id = orderInst.get(pos).spec.id;
                pp.polyPlaced = best.placedGeom;
                pp.angleDeg = best.angleDeg;
                pp.mirrored = false;
                pp.dx = best.dx; pp.dy = best.dy;

                res.placed.add(pp);
                res.placedClear.add(best.clearGeom);

                Envelope ce = best.clearGeom.getEnvelopeInternal();
                res.globMinX = Math.min(res.globMinX, ce.getMinX());
                res.globMinY = Math.min(res.globMinY, ce.getMinY());
                res.globMaxX = Math.max(res.globMaxX, ce.getMaxX());
                res.globMaxY = Math.max(res.globMaxY, ce.getMaxY());

                addLocalAnchors(anchors, ce, minX, minY, maxX, maxY);
            }
        }
        return res;
    }

    //aplicar búsqueda local rápida sobre un individuo usando reinserciones, tweaks de ángulos y swaps.
    static void localSearch(Individual ind,
                        Polygon sheetInset,
                        List<PartInstance> instances,
                        int[][] angleDomain,
                        Random rng,
                        long millisBudget)
    {
        long t0 = System.currentTimeMillis();
        if (ind.layout == null || ind.layout.size() <= 1) return;

        int tries = Math.min(20, ind.layout.size());
        for (int t = 0; t < tries; t++) {
            if (System.currentTimeMillis() - t0 > millisBudget) break;

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

        int swaps = Math.min(25, ind.order.length - 1);
        for (int s = 0; s < swaps; s++) {
            if (System.currentTimeMillis() - t0 > millisBudget) break;
            int i = rng.nextInt(ind.order.length - 1);
            swap(ind.order, i, i + 1);
            evaluate(ind, sheetInset, instances, angleDomain, rng);
        }
    }

    //generar un individuo completamente aleatorio (orden + ángulos) para diversidad en la población.
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

    // crear una semilla “greedy” con orden por área y pequeños jitters para arrancar bien el GA.
    static Individual greedySeed(int N, List<PartInstance> instances, int[][] angleDomain, Random rng){
        Integer[] idx = new Integer[N];
        for (int i=0;i<N;i++) idx[i]=i;

        Arrays.sort(idx, (a,b)-> {
            double aa = instances.get(a).spec.poly.getArea();
            double bb = instances.get(b).spec.poly.getArea();
            double ja = (rng.nextDouble()-0.5)*1e-6*Math.max(1.0, aa);
            double jb = (rng.nextDouble()-0.5)*1e-6*Math.max(1.0, bb);
            return Double.compare((bb+jb), (aa+ja));
        });

        Individual ind = new Individual();
        ind.order = new int[N]; for (int i=0;i<N;i++) ind.order[i] = idx[i];
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

    //seleccionar el mejor individuo de un torneo aleatorio de tamaño k.
    static Individual tournament(List<Individual> pop, Random rng, int k){
        Individual best = null;
        for (int i=0;i<k;i++) {
            Individual c = pop.get(rng.nextInt(pop.size()));
            if (best==null || isBetter(c, best)) best = c;
        }
        return best;
    }

    //cruzar padres con OX para el orden y mezclar ángulos respetando el dominio por gen.
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

    // aplicar mutaciones al orden (swap/reverse) y a los ángulos con tasas adaptadas.
    static void mutate(Individual ch, int[][] angleDomain, Random rng, double mutOrder, double mutAngle){
        if (rng.nextDouble() < mutOrder) {
            int i = rng.nextInt(ch.order.length);
            int j = rng.nextInt(ch.order.length);
            swap(ch.order, i, j);
        }
        if (rng.nextDouble() < mutOrder*0.75) {
            int i = rng.nextInt(ch.order.length-1);
            reverse(ch.order, i, Math.min(ch.order.length-1, i+1+rng.nextInt(Math.max(1, ch.order.length/6))));
        }
        for (int i=0;i<ch.angles.length;i++) {
            if (rng.nextDouble() < mutAngle) {
                int[] dom = angleDomain[ch.order[i]];
                int idx = Arrays.binarySearch(dom, snapToDomain(ch.angles[i], dom));
                if (idx < 0) idx = -idx-1;
                int step = rng.nextBoolean()?1:-1;
                ch.angles[i] = dom[(idx + step + dom.length) % dom.length];
            }
        }
    }

    // evaluar un individuo decodificando su layout y calculando métricas (área y yMax).
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
        ind.area   = totalArea(dr.placed);
        ind.yMax   = usedYmax(dr.placed);
    }

    //comparar dos individuos priorizando mayor área, luego menor yMax y finalmente un jitter estable.
    static boolean isBetter(Individual a, Individual b){
        if ((a.area > b.area)) return true;
        if ((a.area < b.area)) return false;
        if (Math.abs(a.area-b.area) < 1e-9) {
            if (a.yMax < b.yMax) return true;
            if (a.yMax > b.yMax) return false;
            return Long.compare(a.tieNoise, b.tieNoise) < 0;
        }
        return false;
    }

    //  obtener el mejor individuo de una población según el criterio isBetter.
    static Individual bestOf(List<Individual> pop){
        Individual best = pop.get(0);
        for (int i=1;i<pop.size();i++) if (isBetter(pop.get(i), best)) best = pop.get(i);
        return best;
    }

    // generar anclas iniciales (esquinas y grilla) para explorar posiciones de colocación.
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

    // añadir anclas alrededor de la última pieza colocada, filtrando las que salgan del borde.
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

    // muestrear hasta 'max' anclas aleatorias para reducir el costo de evaluación.
    static List<Anchor> sampleAnchors(List<Anchor> anchors, Random rng, int max){
        if (anchors.size() <= max) return new ArrayList<>(anchors);
        ArrayList<Anchor> out = new ArrayList<>(max);
        ArrayList<Integer> idx = new ArrayList<>(anchors.size());
        for (int i=0;i<anchors.size();i++) idx.add(i);
        Collections.shuffle(idx, rng);
        for (int i=0;i<max;i++) out.add(anchors.get(idx.get(i)));
        return out;
    }

    //  detectar solape real (no cuenta contactos tangenciales de borde).
    static boolean collideStrict(Geometry a, Geometry b) {
        return a.intersects(b) && !a.touches(b);
    }

    // micro-compactar una geometría clara moviéndola abajo/izquierda/diagonal si sigue siendo válida.
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

    //  puntuar una colocación favoreciendo menor Y, luego menor X y menor crecimiento del bbox global.
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

    // sumar el área de todas las piezas colocadas.
    static double totalArea(List<PlacedPart> placed){
        double s=0; for (PlacedPart p : placed) s += p.polyPlaced.getArea(); return s;
    }

    // obtener la Y máxima usada por el layout (útil como desempate/compactación vertical).
    static double usedYmax(List<PlacedPart> placed){
        double y=0; for (PlacedPart p : placed) y = Math.max(y, p.polyPlaced.getEnvelopeInternal().getMaxY()); return y;
    }

    // implementar el cruce OX clásico entre dos órdenes.
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

    // máscara para saber de qué padre provino cada posición tras OX (útil para mezclar ángulos).
    static boolean[] originMaskOX(int[] p1, int[] p2, int[] child){
        boolean[] mask = new boolean[child.length];
        for (int i=0;i<child.length;i++) mask[i] = (child[i]==p1[i]);
        return mask;
    }

    // buscar el índice de un valor en un arreglo (lineal, pequeño N).
    static int indexOf(int[] arr, int val){
        for (int i=0;i<arr.length;i++) if (arr[i]==val) return i;
        return -1;
    }

    // verificar si un valor existe en el arreglo.
    static boolean contains(int[] arr, int val){
        for (int v: arr) if (v==val) return true;
        return false;
    }

    // intercambiar dos posiciones en un arreglo.
    static void swap(int[] a, int i, int j){ int t=a[i]; a[i]=a[j]; a[j]=t; }

    // invertir un segmento [i,j] en el arreglo.
    static void reverse(int[] a, int i, int j){
        while (i<j){ swap(a,i++,j--); }
    }

    // barajar un arreglo in-place (Fisher-Yates).
    static void shuffle(int[] a, Random rng){
        for (int i=a.length-1;i>0;i--){
            int j=rng.nextInt(i+1);
            swap(a,i,j);
        }
    }

    // normalizar un ángulo a [0,360).
    static int norm(int ang){ int x = ang % 360; if (x<0) x+=360; return x; }

    // ajustar un ángulo al valor más cercano dentro de un dominio discreto.
    static int snapToDomain(int ang, int[] domain){
        if (domain==null || domain.length==0) return ang;
        int best=domain[0], bestd = Math.abs(norm(ang)-domain[0]);
        for (int v: domain){
            int d = Math.abs(norm(ang)-v);
            if (d<bestd){ best=v; bestd=d; }
        }
        return best;
    }

    //expandir PartSpec en instancias individuales según su cantidad (qty).
    static List<PartInstance> expandInstances(List<PartSpec> specs) {
        List<PartInstance> out = new ArrayList<>();
        for (PartSpec s : specs) {
            int q = Math.max(1, s.qty);
            for (int i=0;i<q;i++) out.add(new PartInstance(s, i));
        }
        return out;
    }

    // construir el dominio de ángulos permitido por cada PartSpec según su step de rotación.
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

    //llenar mediante raster cuando todas las piezas son casi iguales y la plancha es rectangular.
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

                boolean placedHere = false;
                for (PartInstance pi : ord){
                    Integer rem = remaining.get(pi.spec.id);
                    if (rem==null || rem<=0) continue;

                    int angDeg = 0;
                    int idx = instances.indexOf(pi);
                    if (idx >= 0 && idx < angles.length) angDeg = norm(angles[idx]);

                    Geometry g0  = GeomUtils.orient(pi.spec.poly, angDeg, false);
                    Geometry gc0 = (pi.spec.clearance > 0) ? g0.buffer(pi.spec.clearance, 8) : g0;

                    Envelope pe0 = gc0.getEnvelopeInternal();
                    if (Math.abs(pe0.getWidth() - cellW) > cellW*0.03 ||
                        Math.abs(pe0.getHeight()- cellH) > cellH*0.03) continue;

                    double dxa = anchorX - pe0.getMinX();
                    double dya = anchorY - pe0.getMinY();

                    Geometry ccTry = GeomUtils.translate(gc0, dxa, dya);
                    if (!sheetInset.covers(ccTry)) continue;

                    boolean clash = false;
                    for (Geometry pg : placedClear) if (collideStrict(pg, ccTry)) { clash = true; break; }
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
            }
        }
    }

    // calcular coeficiente de variación para decidir si aplica modo rejilla.
    static double cv(double[] a){
        if (a==null || a.length==0) return 1.0;
        double s=0; for(double x:a) s+=x;
        double m=s/a.length; if (m==0) return 1.0;
        double v=0; for(double x:a) v+=(x-m)*(x-m);
        v/=a.length; return Math.sqrt(v)/m;
    }

    // obtener la mediana de un arreglo (copiando y ordenando).
    static double median(double[] a){
        if (a==null || a.length==0) return 0.0;
        double[] cp=a.clone(); Arrays.sort(cp); return cp[cp.length/2];
    }

    // abrir el logger a archivo y arrancar el tee para consola + archivo.
    static void initLogger() {
        try {
            tee = new PrintStream(new FileOutputStream(LOG_FILE, false), true, StandardCharsets.UTF_8);
        } catch (Exception e) {
            tee = null;
        }
        log("[" + TS.format(new Date()) + "] Logger abierto");
    }

    //cerrar el stream de logging si está activo.
    static void closeLogger() { if (tee!=null) tee.close(); }

    //imprimir mensaje con timestamp a consola y archivo (si está habilitado).
    static void log(String s) {
        String msg = TS.format(new Date()) + " | " + s;
        System.out.println(msg);
        if (tee!=null) tee.println(msg);
    }

    // clonar superficialmente un individuo (referencia a layout compartida).
    static Individual cloneShallow(Individual a){
        Individual b = new Individual();
        b.order = a.order.clone();
        b.angles = a.angles.clone();
        b.layout = a.layout;
        b.area = a.area; b.yMax = a.yMax;
        b.tieNoise = a.tieNoise;
        return b;
    }

    //clonar profundamente un individuo (copia de lista de layout).
    static Individual cloneDeep(Individual a){
        Individual b = new Individual();
        b.order = a.order.clone();
        b.angles = a.angles.clone();
        if (a.layout != null) b.layout = new ArrayList<>(a.layout);
        b.area = a.area; b.yMax = a.yMax;
        b.tieNoise = a.tieNoise;
        return b;
    }

    //generar una semilla pseudoaleatoria combinando nanoTime y Math.random.
    static long randomSeed() {
        long s = System.nanoTime() ^ Double.doubleToLongBits(Math.random());
        return (s == 0L ? 1L : s);
    }
}
