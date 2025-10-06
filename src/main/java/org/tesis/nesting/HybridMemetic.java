package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.util.AffineTransformation;
import org.locationtech.jts.operation.union.UnaryUnionOp;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

/**
 * =====================  HybridMemetic  =====================
 * Memético diseñado para integrarse con un algoritmo híbrido.
 * Debido a que MemeticPlacer solo recibe los datos y devuelve un resultado
 * Se programó otro código memético que recibe los datos iniciales para
 * poder usarse en el híbrido. Y recibir como soluciones base las generadas por
 * el GRASP
 */
public class HybridMemetic {

    // ===================== Parámetros GA =====================
    private static final int    POP_SIZE          = 36;
    private static final int    ELITE_SIZE        = 4;
    private static final int    GEN_MAX           = 200;
    private static final int    TOURNAMENT_K      = 3;

    // Mutación 
    private static final double MUT_DX_MAX        = 3.0;  
    private static final double MUT_DY_MAX        = 3.0;
    private static final double MUT_ROT_MAX_DEG   = 2.0;   

    // Búsqueda local
    private static final int    LS_PASSES         = 2;
    private static final double SLIDE_STEP        = 1.0;
    private static final double ROT_STEP_DEG      = 1.0;

    // Formatos y logging
    private static final DecimalFormat DF3  = new DecimalFormat("#,##0.000");
    private static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");


    public static List<PlacedPart> improve(
            Polygon sheet,
            List<PartSpec> parts,
            String logFile,
            long timeLimitMs,
            long seed,
            List<PlacedPart>... seedSolutions) throws Exception {

        List<List<PlacedPart>> seeds = new ArrayList<>();
        if (seedSolutions != null) {
            for (List<PlacedPart> s : seedSolutions) {
                if (s != null) seeds.add(cloneSolution(s));
            }
        }
        return improveFromPopulation(sheet, parts, seeds, logFile, timeLimitMs, seed);
    }

    public static List<PlacedPart> improveFromPopulation(
            Polygon sheet,
            List<PartSpec> parts,
            List<List<PlacedPart>> seedPopulation,
            String logFile,
            long timeLimitMs,
            long seed) throws Exception {

        Objects.requireNonNull(sheet, "sheet no puede ser null");
        Objects.requireNonNull(parts, "parts no puede ser null");

        Random rnd = new Random(seed);
        long t0 = System.currentTimeMillis();

        try (PrintStream log = openLog(logFile)) {
            head(log, "==== HYBRID-MEMETIC | START ====");
            line(log, "CFG | POP_SIZE=" + POP_SIZE + " ELITE=" + ELITE_SIZE + " GEN_MAX=" + GEN_MAX + " K=" + TOURNAMENT_K);
            line(log, "CFG | seed=" + seed + " parts=" + parts.size() + " timeLimitMs=" + timeLimitMs);

            // 1) Población inicial
            List<List<PlacedPart>> population = initPopulation(sheet, parts, seedPopulation, rnd, log);
            List<Individual> pop = evalAll(sheet, population);

            Individual best = bestOf(pop);

            // 2) Evolución
            int gen = 0;
            while (gen < GEN_MAX && !timeout(t0, timeLimitMs)) {
                List<Individual> next = new ArrayList<>(POP_SIZE);

                // Elitismo
                next.addAll(selectElites(pop, ELITE_SIZE));

                // Reproducción
                while (next.size() < POP_SIZE) {
                    Individual p1 = tournament(pop, rnd);
                    Individual p2 = tournament(pop, rnd);
                    List<PlacedPart> childSol = crossover(sheet, p1.solution, p2.solution, rnd);
                    mutate(sheet, childSol, rnd);
                    localSearch(sheet, childSol);
                    next.add(eval(sheet, childSol));
                }

                pop = next;
                Individual genBest = bestOf(pop);
                if (compare(genBest, best) < 0) best = genBest;

                if ((gen % 10 == 0) || timeout(t0, timeLimitMs) || gen == GEN_MAX - 1) {
                    line(log, "GEN " + gen + " | bestPlaced=" + best.placedCount +
                            " | area=" + DF3.format(best.areaPlaced) +
                            " | yMax=" + DF3.format(best.yMax));
                }
                gen++;
            }

            // 3) Resumen
            long t1 = System.currentTimeMillis();
            long totalMs = t1 - t0;

            summary(log, "RESUMEN (HYBRID-MEMETIC)");
            line(log, "Seed               : " + seed);
            line(log, "Piezas a posicionar: " + parts.size());
            line(log, "Piezas colocadas   : " + best.placedCount);
            line(log, "% aprovechamiento  : " + DF3.format(best.utilizationPct) + " %");
            line(log, "Área colocada      : " + DF3.format(best.areaPlaced));
            line(log, "Área plancha       : " + DF3.format(best.areaSheet));
            line(log, "Área desperdiciada : " + DF3.format(best.areaWaste));
            line(log, "Ymax               : " + DF3.format(best.yMax));
            line(log, "Tiempo total       : " + totalMs + " ms (" + DF3.format(totalMs / 1000.0) + " s)");
            dash(log);
            line(log, "BEST | HYB-MEM | colocadas=" + best.placedCount +
                    " | area=" + DF3.format(best.areaPlaced) +
                    " | yMax=" + DF3.format(best.yMax));
            dash(log);

            return best.solution;
        }
    }

    // ===================== Inicialización =====================

    private static List<List<PlacedPart>> initPopulation(
            Polygon sheet,
            List<PartSpec> parts,
            List<List<PlacedPart>> seeds,
            Random rnd,
            PrintStream log) {

        List<List<PlacedPart>> pop = new ArrayList<>(POP_SIZE);

        // 1) Inyectar seeds (clonadas y saneadas)
        if (seeds != null) {
            for (List<PlacedPart> s : seeds) {
                if (s == null) continue;
                List<PlacedPart> clean = sanitizeSolution(sheet, s);
                if (!clean.isEmpty()) pop.add(clean);
            }
        }

        // 2) Completar con variantes "jitter" a partir de lo ya inyectado
        while (pop.size() < POP_SIZE) {
            List<PlacedPart> base = (!pop.isEmpty()) ? pop.get(rnd.nextInt(pop.size())) : new ArrayList<>();
            List<PlacedPart> variant = jitterVariant(sheet, base, rnd);
            localSearch(sheet, variant);
            pop.add(variant);
        }

        line(log, "INIT | población=" + pop.size() + " (seeds=" + (seeds != null ? seeds.size() : 0) + ")");
        return pop;
    }

    private static List<PlacedPart> sanitizeSolution(Polygon sheet, List<PlacedPart> sol) {
        List<PlacedPart> out = new ArrayList<>();
        Geometry acc = null;
        for (PlacedPart p : sol) {
            if (p == null || p.polyPlaced == null || p.polyPlaced.isEmpty()) continue;
            if (!p.polyPlaced.within(sheet)) continue;
            if (acc == null) {
                acc = p.polyPlaced;
                out.add(clonePlaced(p));
            } else if (!p.polyPlaced.intersects(acc)) {
                out.add(clonePlaced(p));
                acc = UnaryUnionOp.union(Arrays.asList(acc, p.polyPlaced));
            }
        }
        return out;
    }

    private static List<PlacedPart> jitterVariant(Polygon sheet, List<PlacedPart> base, Random rnd) {
        List<PlacedPart> out = new ArrayList<>(base.size());
        Geometry acc = null;
        for (PlacedPart p : base) {
            PlacedPart q = clonePlaced(p);
            for (int tries = 0; tries < 8; tries++) {
                double dx = (rnd.nextDouble() * 2 - 1) * MUT_DX_MAX;
                double dy = (rnd.nextDouble() * 2 - 1) * MUT_DY_MAX;
                double drot = (rnd.nextDouble() * 2 - 1) * MUT_ROT_MAX_DEG;
                Geometry moved = affineTransform(q.polyPlaced, dx, dy, Math.toRadians(drot));
                if (moved == null || moved.isEmpty()) continue;
                if (!moved.within(sheet)) continue;
                if (acc != null && moved.intersects(acc)) continue;
                q.polyPlaced = moved;
                q.dx += dx; q.dy += dy; q.angleDeg += drot;
                break;
            }
            if (q.polyPlaced != null && !q.polyPlaced.isEmpty() && q.polyPlaced.within(sheet)) {
                if (acc == null) {
                    acc = q.polyPlaced;
                    out.add(q);
                } else if (!q.polyPlaced.intersects(acc)) {
                    out.add(q);
                    acc = UnaryUnionOp.union(Arrays.asList(acc, q.polyPlaced));
                }
            }
        }
        return out;
    }

    // ===================== Evaluación / comparación =====================

    private static Individual eval(Polygon sheet, List<PlacedPart> sol) {
        Individual ind = new Individual();
        ind.solution = sanitizeSolution(sheet, sol);
        ind.placedCount = ind.solution.size();
        ind.areaSheet = sheet.getArea();
        ind.areaPlaced = placedArea(ind.solution);
        ind.areaWaste  = Math.max(0.0, ind.areaSheet - ind.areaPlaced);
        ind.utilizationPct = ind.areaSheet <= 0 ? 0.0 : (ind.areaPlaced / ind.areaSheet) * 100.0;
        ind.yMax = yMax(ind.solution);
        ind.elapsedMs = 0;
        return ind;
    }

    private static List<Individual> evalAll(Polygon sheet, List<List<PlacedPart>> sols) {
        List<Individual> out = new ArrayList<>(sols.size());
        for (List<PlacedPart> s : sols) out.add(eval(sheet, s));
        return out;
    }

    /** -1 si a es mejor que b; +1 si peor; 0 empate. */
    private static int compare(Individual a, Individual b) {
        if (a == null && b == null) return 0;
        if (a == null) return 1;
        if (b == null) return -1;

        if (a.placedCount != b.placedCount) return (a.placedCount > b.placedCount) ? -1 : 1;
        if (Double.compare(a.areaPlaced, b.areaPlaced) != 0) return (a.areaPlaced > b.areaPlaced) ? -1 : 1;
        if (Double.compare(a.yMax, b.yMax) != 0) return (a.yMax < b.yMax) ? -1 : 1;
        return 0;
    }

    private static Individual bestOf(List<Individual> pop) {
        Individual best = null;
        for (Individual i : pop) if (compare(i, best) < 0) best = i;
        return best;
    }

    private static List<Individual> selectElites(List<Individual> pop, int k) {
        List<Individual> copy = new ArrayList<>(pop);
        copy.sort(HybridMemetic::compare);
        return new ArrayList<>(copy.subList(0, Math.min(k, copy.size())));
    }

    private static Individual tournament(List<Individual> pop, Random rnd) {
        Individual best = null;
        for (int i = 0; i < TOURNAMENT_K; i++) {
            Individual cand = pop.get(rnd.nextInt(pop.size()));
            if (compare(cand, best) < 0) best = cand;
        }
        return best;
    }

    private static double placedArea(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        List<Geometry> gs = new ArrayList<>(placed.size());
        for (PlacedPart p : placed) {
            if (p != null && p.polyPlaced != null && !p.polyPlaced.isEmpty()) gs.add(p.polyPlaced);
        }
        if (gs.isEmpty()) return 0.0;
        Geometry u = UnaryUnionOp.union(gs);
        return (u != null) ? u.getArea() : 0.0;
    }

    private static double yMax(List<PlacedPart> placed) {
        if (placed == null || placed.isEmpty()) return 0.0;
        double ymax = Double.NEGATIVE_INFINITY;
        for (PlacedPart p : placed) {
            if (p == null || p.polyPlaced == null || p.polyPlaced.isEmpty()) continue;
            Envelope env = p.polyPlaced.getEnvelopeInternal();
            if (env != null) ymax = Math.max(ymax, env.getMaxY());
        }
        return (ymax == Double.NEGATIVE_INFINITY) ? 0.0 : ymax;
    }

    // ===================== Operadores genéticos =====================

    /** Crossover por pieza: toma la mejor variante válida de cada padre, priorizando menor Ymax incremental. */
    private static List<PlacedPart> crossover(
            Polygon sheet,
            List<PlacedPart> a,
            List<PlacedPart> b,
            Random rnd) {

        Map<String, List<PlacedPart>> cand = new HashMap<>();
        for (PlacedPart p : a) addCand(cand, p);
        for (PlacedPart p : b) addCand(cand, p);

        List<String> ids = new ArrayList<>(cand.keySet());
        ids.sort((i1, i2) -> Integer.compare(cand.get(i2).size(), cand.get(i1).size())); // primero las que tienen más opciones

        List<PlacedPart> child = new ArrayList<>();
        Geometry acc = null;

        for (String id : ids) {
            List<PlacedPart> opts = cand.get(id);
            Collections.shuffle(opts, rnd);

            PlacedPart best = null;
            double bestScore = Double.POSITIVE_INFINITY;

            for (PlacedPart opt : opts) {
                if (opt == null || opt.polyPlaced == null || opt.polyPlaced.isEmpty()) continue;
                if (!opt.polyPlaced.within(sheet)) continue;
                if (acc != null && opt.polyPlaced.intersects(acc)) continue;

                double yIf = yMaxIfAdded(child, opt);
                if (yIf < bestScore) {
                    bestScore = yIf;
                    best = opt;
                }
            }

            if (best != null) {
                PlacedPart q = clonePlaced(best);
                child.add(q);
                acc = (acc == null) ? q.polyPlaced : UnaryUnionOp.union(Arrays.asList(acc, q.polyPlaced));
            }
        }
        return child;
    }

    private static void mutate(Polygon sheet, List<PlacedPart> sol, Random rnd) {
        if (sol.isEmpty()) return;
        int toMut = Math.max(1, sol.size() / 6);
        Collections.shuffle(sol, rnd);

        for (int i = 0; i < toMut; i++) {
            PlacedPart p = sol.get(i);
            Geometry accWithout = unionExcept(p, sol);
            for (int tries = 0; tries < 10; tries++) {
                double dx = (rnd.nextDouble() * 2 - 1) * MUT_DX_MAX;
                double dy = (rnd.nextDouble() * 2 - 1) * MUT_DY_MAX;
                double drot = (rnd.nextDouble() * 2 - 1) * MUT_ROT_MAX_DEG;
                Geometry moved = affineTransform(p.polyPlaced, dx, dy, Math.toRadians(drot));
                if (moved == null || moved.isEmpty()) continue;
                if (!moved.within(sheet)) continue;
                if (accWithout != null && moved.intersects(accWithout)) continue;

                p.polyPlaced = moved;
                p.dx += dx; p.dy += dy; p.angleDeg += drot;
                break;
            }
        }
    }

    // ===================== Búsqueda local =====================

    private static void localSearch(Polygon sheet, List<PlacedPart> sol) {
        if (sol.isEmpty()) return;
        for (int pass = 0; pass < LS_PASSES; pass++) {
            sol.sort(Comparator.comparingDouble(HybridMemetic::minYOf));
            for (PlacedPart p : sol) {
                slidePiece(sheet, p, sol, 0, -SLIDE_STEP); // hacia abajo
                slidePiece(sheet, p, sol, -SLIDE_STEP, 0); // hacia la izquierda
                microRotateImprove(sheet, p, sol);
            }
        }
    }

    private static void slidePiece(Polygon sheet, PlacedPart p, List<PlacedPart> sol, double stepx, double stepy) {
        if (p.polyPlaced == null || p.polyPlaced.isEmpty()) return;
        for (int k = 0; k < 50; k++) {
            Geometry candidate = affineTransform(p.polyPlaced, stepx, stepy, 0);
            if (candidate == null || candidate.isEmpty()) break;
            if (!candidate.within(sheet)) break;
            if (intersectsOthers(candidate, p, sol)) break;
            p.polyPlaced = candidate;
            p.dx += stepx;
            p.dy += stepy;
        }
    }

    private static void microRotateImprove(Polygon sheet, PlacedPart p, List<PlacedPart> sol) {
        if (p.polyPlaced == null || p.polyPlaced.isEmpty()) return;
        double currentY = yMax(sol);
        double[] deltas = new double[]{Math.toRadians(ROT_STEP_DEG), Math.toRadians(-ROT_STEP_DEG)};
        for (double drot : deltas) {
            Geometry candidate = affineTransform(p.polyPlaced, 0, 0, drot);
            if (candidate == null || candidate.isEmpty()) continue;
            if (!candidate.within(sheet)) continue;
            if (intersectsOthers(candidate, p, sol)) continue;
            Geometry old = p.polyPlaced;
            p.polyPlaced = candidate;
            double newY = yMax(sol);
            if (newY <= currentY) {
                p.angleDeg += Math.toDegrees(drot);
                currentY = newY;
            } else {
                p.polyPlaced = old;
            }
        }
    }

    // ===================== Utilidades geométricas =====================

    /** Rotación sobre el centroide + traslación. */
    private static Geometry affineTransform(Geometry g, double dx, double dy, double rotRad) {
        if (g == null || g.isEmpty()) return g;
        Coordinate c = g.getCentroid().getCoordinate();
        AffineTransformation at = new AffineTransformation();
        if (rotRad != 0) at.rotate(rotRad, c.x, c.y);
        if (dx != 0 || dy != 0) at.translate(dx, dy);
        return at.transform(g);
    }

    private static Geometry unionExcept(PlacedPart except, List<PlacedPart> sol) {
        List<Geometry> gs = new ArrayList<>();
        for (PlacedPart q : sol) {
            if (q == null || q.polyPlaced == null || q.polyPlaced.isEmpty()) continue;
            if (q == except) continue;
            gs.add(q.polyPlaced);
        }
        if (gs.isEmpty()) return null;
        return UnaryUnionOp.union(gs);
    }

    private static boolean intersectsOthers(Geometry geom, PlacedPart self, List<PlacedPart> sol) {
        for (PlacedPart q : sol) {
            if (q == null || q == self || q.polyPlaced == null || q.polyPlaced.isEmpty()) continue;
            if (geom.intersects(q.polyPlaced)) return true;
        }
        return false;
    }

    private static double yMaxIfAdded(List<PlacedPart> base, PlacedPart toAdd) {
        double y = yMax(base);
        if (toAdd != null && toAdd.polyPlaced != null && !toAdd.polyPlaced.isEmpty()) {
            Envelope e = toAdd.polyPlaced.getEnvelopeInternal();
            if (e != null) y = Math.max(y, e.getMaxY());
        }
        return y;
    }

    private static double minYOf(PlacedPart p) {
        if (p == null || p.polyPlaced == null || p.polyPlaced.isEmpty()) return 0.0;
        Envelope e = p.polyPlaced.getEnvelopeInternal();
        return (e == null) ? 0.0 : e.getMinY();
    }

    private static void addCand(Map<String, List<PlacedPart>> cand, PlacedPart p) {
        if (p == null || p.id == null) return;
        cand.computeIfAbsent(p.id, k -> new ArrayList<>()).add(clonePlaced(p));
    }

    private static PlacedPart clonePlaced(PlacedPart p) {
        PlacedPart q = new PlacedPart();
        q.id = p.id;
        q.polyPlaced = (p.polyPlaced != null) ? (Geometry) p.polyPlaced.copy() : null;
        q.angleDeg = p.angleDeg;
        q.mirrored = p.mirrored;
        q.dx = p.dx;
        q.dy = p.dy;
        return q;
    }

    private static List<PlacedPart> cloneSolution(List<PlacedPart> s) {
        List<PlacedPart> out = new ArrayList<>(s.size());
        for (PlacedPart p : s) out.add(clonePlaced(p));
        return out;
    }

    private static boolean timeout(long t0, long limitMs) {
        if (limitMs <= 0) return false;
        return System.currentTimeMillis() - t0 >= limitMs;
    }

    // ===================== Logging =====================

    private static PrintStream openLog(String logFile) throws IOException {
        if (logFile == null || logFile.isEmpty()) {
            new File("out").mkdirs();
            logFile = "out/LogHybrid.log";
        } else {
            File lf = new File(logFile).getAbsoluteFile();
            File parent = lf.getParentFile();
            if (parent != null) parent.mkdirs();
        }
        // PrintStream(OutputStream out, boolean autoFlush, String encoding)
        return new PrintStream(new FileOutputStream(logFile, true), true, StandardCharsets.UTF_8.name());
    }

    private static void head(PrintStream log, String title) {
        dash(log);
        line(log, title);
    }

    private static void summary(PrintStream log, String title) {
        dash(log);
        line(log, title);
    }

    private static void dash(PrintStream log) {
        line(log, "------------------------------");
    }

    private static void line(PrintStream log, String msg) {
        log.println(TS.format(new Date()) + " | " + msg);
    }

    // ===================== POJO Individual =====================

    private static class Individual {
        List<PlacedPart> solution;
        int    placedCount;
        double areaPlaced;
        double areaSheet;
        double areaWaste;
        double yMax;
        double utilizationPct;
        long   elapsedMs;
    }
}
