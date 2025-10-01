package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.util.*;

public class BottomLeftPlacer {

    // función que coloca las piezas usando el barrido bottom-left
    static List<PlacedPart> place(Polygon sheetInset, List<PartSpec> specs) {
        List<PlacedPart> placed = new ArrayList<>();
        List<Geometry> placedClearance = new ArrayList<>();

        Envelope env = sheetInset.getEnvelopeInternal();
        double width  = env.getWidth();
        double height = env.getHeight();
        double gridStep = Math.max(1.0, Math.min(width, height) * 0.01); // paso de la grilla

        for (PartSpec sp : specs) {
            int rotStep = Math.max(1, sp.rotationStepDeg);

            List<Integer> angles = new ArrayList<>();
            for (int a = 0; a < 360; a += rotStep) angles.add(a);

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

                                PlacedPart pp = new PlacedPart();
                                pp.id = sp.id;
                                pp.polyPlaced = candidate;
                                pp.angleDeg = aDeg;
                                pp.mirrored = mir;
                                pp.dx = dx;
                                pp.dy = dy;

                                placed.add(pp);
                                placedClearance.add(candidateClear);
                                ok = true;
                                break outer;
                            }
                        }
                    }
                }
                if (!ok) {
                    System.err.println("Aviso: no se pudo colocar una copia de " + sp.id + " (" + (copy + 1) + "/" + sp.qty + ")");
                }
            }
        }
        return placed;
    }
    //
}
