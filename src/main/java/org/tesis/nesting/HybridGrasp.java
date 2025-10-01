package org.tesis.nesting;

import org.locationtech.jts.geom.Polygon;

import java.io.PrintStream;
import java.text.DecimalFormat;
import java.text.SimpleDateFormat;
import java.util.*;

public final class HybridGrasp {

    private static final SimpleDateFormat TS = new SimpleDateFormat("HH:mm:ss.SSS");
    private static final DecimalFormat DF = new DecimalFormat("#,##0.###");

    private HybridGrasp() {}

    // corre GRASP varias veces con diferentes seeds y devuelve la lista de layouts
    public static List<List<PlacedPart>> multiStart(Polygon sheetInset,
                                                    List<PartSpec> specs,
                                                    int k,
                                                    long baseSeed,
                                                    PrintStream log) {
        Objects.requireNonNull(sheetInset, "sheetInset");
        Objects.requireNonNull(specs, "specs");
        List<List<PlacedPart>> seeds = new ArrayList<>(k);
        for (int i = 0; i < k; i++) {
            long seed = baseSeed + i;
            long t0 = System.currentTimeMillis();
            List<PlacedPart> layout = callGrasp(sheetInset, specs, seed, log);
            int placed = (layout == null ? 0 : layout.size());
            if (log != null) log.println(tag() + String.format("GRASP seed %02d (seed=%d) -> placed=%d (%.3fs)",
                    i, seed, placed, (System.currentTimeMillis() - t0) / 1000.0));
            seeds.add(layout == null ? Collections.emptyList() : layout);
        }
        return seeds;
    }

    // elige el mejor layout entre los candidatos (más piezas colocadas o por scorer)
    public static List<PlacedPart> pickBest(List<List<PlacedPart>> candidates,
                                            LayoutScorer scorer,
                                            PrintStream log) {
        List<PlacedPart> best = Collections.emptyList();
        double bestScore = Double.NEGATIVE_INFINITY;
        int bestPlaced = -1;

        for (int i = 0; i < candidates.size(); i++) {
            List<PlacedPart> lay = (candidates.get(i) == null ? Collections.emptyList() : candidates.get(i));
            int placed = lay.size();
            double s = (scorer != null ? scorer.score(lay) : placed);
            boolean better = (placed > bestPlaced) || (placed == bestPlaced && s > bestScore);
            if (better) {
                best = lay; bestPlaced = placed; bestScore = s;
            }
            if (log != null) log.println(tag() + String.format("CAND %02d -> placed=%d score=%s",
                    i, placed, DF.format(s)));
        }
        if (log != null) log.println(tag() + String.format("BEST_GRASP -> placed=%d score=%s",
                (best == null ? 0 : best.size()), DF.format(bestScore)));
        return best;
    }

    // usa GRASP para reinsertar piezas en un espacio libre (útil para repair en algoritmos evolutivos)
    public static List<PlacedPart> reinsertOnFree(Polygon freeSpace,
                                                  List<PartSpec> subset,
                                                  long seed,
                                                  PrintStream log) {
        if (freeSpace == null || subset == null || subset.isEmpty()) {
            return Collections.emptyList();
        }
        long t0 = System.currentTimeMillis();
        List<PlacedPart> res = callGrasp(freeSpace, subset, seed, log);
        if (log != null) {
            int placed = (res == null ? 0 : res.size());
            log.println(tag() + String.format("REPAIR_GRASP seed=%d subset=%d -> placed=%d (%.3fs)",
                    seed, subset.size(), placed, (System.currentTimeMillis() - t0) / 1000.0));
        }
        return (res == null ? Collections.emptyList() : res);
    }

    // interfaz para definir un scorer personalizado de layouts
    public interface LayoutScorer {
        double score(List<PlacedPart> layout);
    }

    // invoca métodos de GraspPlacer usando reflexión (soporta varias firmas)
    @SuppressWarnings("unchecked")
    private static List<PlacedPart> callGrasp(Polygon sheetInset,
                                              List<PartSpec> specs,
                                              long seed,
                                              PrintStream log) {
        try {
            return (List<PlacedPart>) GraspPlacer.class
                    .getMethod("place", Polygon.class, List.class, long.class, PrintStream.class)
                    .invoke(null, sheetInset, specs, seed, log);
        } catch (NoSuchMethodException ignore) {
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return Collections.emptyList();
        }

        try {
            return (List<PlacedPart>) GraspPlacer.class
                    .getMethod("place", Polygon.class, List.class, long.class)
                    .invoke(null, sheetInset, specs, seed);
        } catch (NoSuchMethodException ignore) {
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return Collections.emptyList();
        }

        try {
            return (List<PlacedPart>) GraspPlacer.class
                    .getMethod("place", Polygon.class, List.class)
                    .invoke(null, sheetInset, specs);
        } catch (Exception e) {
            if (log != null) e.printStackTrace(log);
            return Collections.emptyList();
        }
    }

    // genera el prefijo de log con timestamp
    private static String tag() {
        return "[" + TS.format(new Date()) + "] ";
    }
}
