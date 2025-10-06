package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.io.File;
import java.io.FileOutputStream;
import java.io.PrintStream;
import java.nio.charset.StandardCharsets;
import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.time.LocalTime;
import java.time.format.DateTimeFormatter;
import java.util.*;

public class BottomLeftPlacer {

    /** Logger con timestamp tipo [HH:mm:ss.SSS] */
    static final class BLLogger implements AutoCloseable {
        private final PrintStream out;
        private final DateTimeFormatter fmt = DateTimeFormatter.ofPattern("HH:mm:ss.SSS");

        BLLogger(PrintStream out) {
            this.out = out;
        }

        static BLLogger toFile(String path) {
            try {
                File f = new File(path);
                File dir = f.getParentFile();
                if (dir != null) dir.mkdirs();
                PrintStream ps = new PrintStream(new FileOutputStream(f, /*append*/false), true, StandardCharsets.UTF_8);
                return new BLLogger(ps);
            } catch (Exception e) {
                throw new RuntimeException("No se pudo abrir el log " + path, e);
            }
        }

        void log(String msg) {
            String t = "[" + LocalTime.now().format(fmt) + "] ";
            out.println(t + msg);
        }

        void logf(String pattern, Object... args) {
            log(String.format(Locale.US, pattern, args));
        }

        @Override public void close() {
            out.flush();
            out.close();
        }
    }

    /* Resultado mínimo para el resumen. */
    public static final class BLReport {
        long seed;
        int piezasAColocar;
        int piezasColocadas;
        double areaColocada;
        double areaPlancha;
        double areaDesperdiciada;
        double yMax;
        long tiempoMs;

        String formatResumen() {
            // Números con miles y 3 decimales, estilo 1,468,792.762
            DecimalFormatSymbols sym = new DecimalFormatSymbols(Locale.US);
            DecimalFormat f3 = new DecimalFormat("#,##0.000", sym);
            DecimalFormat f0 = new DecimalFormat("#,##0", sym);

            String tiempoMsFmt = f0.format(tiempoMs);
            String tiempoSegFmt = f3.format(tiempoMs / 1000.0);
            String utilPct = (areaPlancha > 0) ? f3.format((areaColocada / areaPlancha) * 100.0) : "0.000";

            // Cada línea con su propio timestamp:
            DateTimeFormatter fmt = DateTimeFormatter.ofPattern("HH:mm:ss.SSS");
            String now = "[" + LocalTime.now().format(fmt) + "] ";
            StringBuilder sb = new StringBuilder();
            sb.append(now).append("------------------------------\n");
            sb.append(now).append(now).append("RESUMEN FINAL (HYBRID)\n");
            sb.append(now).append(now).append("Seed               : ").append(seed).append("\n");
            sb.append(now).append(now).append("Piezas a posicionar: ").append(piezasAColocar).append("\n");
            sb.append(now).append(now).append("Piezas colocadas   : ").append(piezasColocadas).append("\n");
            sb.append(now).append(now).append("% aprovechamiento  : ").append(utilPct).append(" %\n");
            sb.append(now).append(now).append("Área colocada      : ").append(f3.format(areaColocada)).append("\n");
            sb.append(now).append(now).append("Área plancha       : ").append(f3.format(areaPlancha)).append("\n");
            sb.append(now).append(now).append("Área desperdiciada : ").append(f3.format(areaDesperdiciada)).append("\n");
            sb.append(now).append(now).append("Ymax               : ").append(f3.format(yMax)).append("\n");
            sb.append(now).append(now).append("Tiempo total       : ").append(tiempoMsFmt).append(" ms (").append(tiempoSegFmt).append(" s)\n");
            sb.append(now).append("------------------------------");
            return sb.toString();
        }
    }

    // función que coloca las piezas usando el barrido bottom-left
    public static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        // Por compatibilidad: crea log por defecto
        try (BLLogger log = BLLogger.toFile("LogBottomLeft.log")) {
            return placeWithLog(sheetInset, specs, log, System.nanoTime());
        }
    }

    /* Variante con logger y semilla explícita (útil para reproducibilidad). */
    public static List<PlacedPart> placeWithLog(Polygon sheetInset, List<PartSpec> specs, BLLogger log, long seed) {
        long t0 = System.currentTimeMillis();
        Random rnd = new Random(seed);

        List<PlacedPart> placed = new ArrayList<>();
        List<Geometry> placedClearance = new ArrayList<>();

        Envelope env = sheetInset.getEnvelopeInternal();
        double width  = env.getWidth();
        double height = env.getHeight();
        double gridStep = Math.max(1.0, Math.min(width, height) * 0.01); // paso de la grilla

        // Métricas en vivo
        final int totalAColocar = specs.stream().mapToInt(sp -> sp.qty).sum();
        double areaColocada = 0.0;
        double areaPlancha  = sheetInset.getArea();
        double yMax = 0.0;

        log.logf("==== Bottom-Left | START ====");
        log.logf("Seed=%d  gridStep=%.4f  (minDim=%.2f)", seed, gridStep, Math.min(width, height));
        log.logf("Piezas a posicionar: %d", totalAColocar);

        int colocadas = 0;
        int intento = 0;

        for (PartSpec sp : specs) {
            int rotStep = Math.max(1, sp.rotationStepDeg);

            List<Integer> angles = new ArrayList<>();
            for (int a = 0; a < 360; a += rotStep) angles.add(a);

            // Pequeño shuffle para diversidad
            Collections.shuffle(angles, rnd);

            List<Boolean> mirrors = sp.mirrorAllowed ? Arrays.asList(false, true) : Collections.singletonList(false);

            for (int copy = 0; copy < sp.qty; copy++) {
                boolean ok = false;

                outer:
                for (boolean mir : mirrors) {
                    for (int aDeg : angles) {
                        Geometry oriented = GeomUtils.orient(sp.poly, aDeg, mir);
                        Geometry orientedClear = (sp.clearance > 0) ? oriented.buffer(sp.clearance, 8) : oriented;

                        Envelope pe = orientedClear.getEnvelopeInternal();
                        double minX = env.getMinX();
                        double minY = env.getMinY();
                        double maxX = env.getMaxX() - pe.getWidth();
                        double maxY = env.getMaxY() - pe.getHeight();
                        if (maxX < minX || maxY < minY) continue;

                        for (double y = minY; y <= maxY + 1e-9; y += gridStep) {
                            for (double x = minX; x <= maxX + 1e-9; x += gridStep) {
                                double dx = x - pe.getMinX();
                                double dy = y - pe.getMinY();

                                Geometry candidateClear = GeomUtils.translate(orientedClear, dx, dy);
                                if (!sheetInset.contains(candidateClear)) continue;

                                boolean intersects = false;
                                for (Geometry gPlaced : placedClearance) {
                                    if (gPlaced.intersects(candidateClear)) {
                                        intersects = true;
                                        break;
                                    }
                                }
                                if (intersects) continue;

                                Geometry candidate = GeomUtils.translate(oriented, dx, dy);

                               // Aceptar
                                PlacedPart pp = new PlacedPart();
                                pp.id        = sp.id;
                                pp.polyPlaced= candidate;   // <-- tu campo real
                                pp.angleDeg  = aDeg;        // double en tu clase, perfecto
                                pp.mirrored  = mir;
                                pp.dx        = dx;          // guarda desplazamiento aplicado
                                pp.dy        = dy;

                                placed.add(pp);
                                placedClearance.add(candidateClear);

                                areaColocada += candidate.getArea();
                                Envelope ce = candidate.getEnvelopeInternal();
                                yMax = Math.max(yMax, ce.getMaxY() - env.getMinY());

                                colocadas++;
                                ok = true;


                                if (colocadas % 5 == 0 || colocadas == totalAColocar) {
                                    double util = (areaPlancha > 0) ? (areaColocada / areaPlancha) * 100.0 : 0.0;
                                    log.logf("Progreso: colocadas=%d/%d  util=%.3f%%  yMax=%.3f",
                                            colocadas, totalAColocar, util, yMax);
                                }
                                break outer;
                            }
                        }
                    }
                }
                if (!ok) {
                    log.logf("Aviso: no se pudo colocar una copia de %s (%d/%d)", sp.id, (copy + 1), sp.qty);
                }
                intento++;
            }
        }

        long t1 = System.currentTimeMillis();

        // Reporte final
        BLReport rep = new BLReport();
        rep.seed = seed;
        rep.piezasAColocar = totalAColocar;
        rep.piezasColocadas = colocadas;
        rep.areaColocada = areaColocada;
        rep.areaPlancha = areaPlancha;
        rep.areaDesperdiciada = Math.max(0.0, areaPlancha - areaColocada);
        rep.yMax = yMax;
        rep.tiempoMs = (t1 - t0);

        log.log(rep.formatResumen());

        return placed;
    }
}
