package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import org.locationtech.jts.geom.util.AffineTransformation;

import java.util.*;
import java.util.stream.Collectors;

public class GeomUtils {

    // construye el polígono de la plancha a partir de las filas con entidad SHEET
    static Polygon buildSheet(List<VertexRow> rows, GeometryFactory gf) {
        Map<String, List<VertexRow>> g = rows.stream()
                .filter(r -> "SHEET".equalsIgnoreCase(r.entityType) && r.ring == 0)
                .collect(Collectors.groupingBy(r -> r.entityId));
        if (g.isEmpty()) return null;
        String sheetId = g.keySet().iterator().next();
        List<VertexRow> pts = g.get(sheetId);
        pts.sort(Comparator.comparingInt(r -> r.pointIdx));
        LinearRing shell = gf.createLinearRing(toClosedCoords(pts));
        return gf.createPolygon(shell);
    }

    // construye la lista de piezas (PartSpec) a partir de las filas con entidad PART
    static List<PartSpec> buildParts(List<VertexRow> rows, GeometryFactory gf) {
        Map<String, List<VertexRow>> outer = rows.stream()
                .filter(r -> "PART".equalsIgnoreCase(r.entityType) && r.ring == 0)
                .collect(Collectors.groupingBy(r -> r.entityId));

        List<PartSpec> list = new ArrayList<>();
        for (Map.Entry<String, List<VertexRow>> e : outer.entrySet()) {
            String id = e.getKey();
            List<VertexRow> pts = e.getValue().stream()
                    .sorted(Comparator.comparingInt(r -> r.pointIdx))
                    .toList();

            LinearRing shell = gf.createLinearRing(toClosedCoords(pts));

            // si la pieza tiene huecos, agruparlos por ring y crearlos como LinearRing[]
            Map<Integer, List<VertexRow>> byRing = rows.stream()
                    .filter(r -> "PART".equalsIgnoreCase(r.entityType) && r.entityId.equals(id) && r.ring > 0)
                    .collect(Collectors.groupingBy(r -> r.ring));
            LinearRing[] holes = null;
            if (!byRing.isEmpty()) {
                holes = new LinearRing[byRing.size()];
                int i = 0;
                for (int ring : new TreeSet<>(byRing.keySet())) {
                    List<VertexRow> hp = byRing.get(ring).stream()
                            .sorted(Comparator.comparingInt(r -> r.pointIdx)).toList();
                    holes[i++] = gf.createLinearRing(toClosedCoords(hp));
                }
            }

            Polygon poly = gf.createPolygon(shell, holes);

            PartSpec sp = new PartSpec();
            sp.id = id;
            sp.poly = poly;

            Integer qty = firstNonNull(pts, r -> r.qty);
            Integer rot = firstNonNull(pts, r -> r.rotationStepDeg);
            Boolean mir = firstNonNull(pts, r -> r.mirrorAllowed);
            Double  clr = firstNonNull(pts, r -> r.clearance);
            Double  gra = firstNonNull(pts, r -> r.grainAngleDeg);

            if (qty != null && qty > 0) sp.qty = qty;
            if (rot != null && rot > 0) sp.rotationStepDeg = rot;
            if (mir != null) sp.mirrorAllowed = mir;
            if (clr != null && clr >= 0) sp.clearance = clr;
            if (gra != null) sp.grainAngleDeg = gra;

            list.add(sp);
        }
        return list;
    }

    // convierte la lista de vértices en un arreglo de coordenadas cerrado (último = primero)
    static Coordinate[] toClosedCoords(List<VertexRow> pts) {
        if (pts.size() < 3) throw new IllegalArgumentException("Polígono requiere >=3 puntos");
        Coordinate[] c = new Coordinate[pts.size() + 1];
        for (int i=0;i<pts.size();i++) c[i] = new Coordinate(pts.get(i).x, pts.get(i).y);
        c[c.length-1] = new Coordinate(pts.get(0).x, pts.get(0).y);
        return c;
    }

    // devuelve el primer valor no nulo encontrado en la lista aplicando un getter
    static <T> T firstNonNull(List<VertexRow> list, java.util.function.Function<VertexRow, T> f) {
        for (VertexRow r : list) { T v = f.apply(r); if (v != null) return v; }
        return null;
    }

    // genera una plancha insetada (offset negativo). Si genera varios polígonos, escoge el mayor.
    static Polygon insetSheet(Polygon sheet, double inset) {
        if (inset <= 0) return sheet;
        Geometry g = sheet.buffer(-inset, 8);
        if (g instanceof Polygon) return (Polygon) g;
        if (g instanceof MultiPolygon mp && mp.getNumGeometries() > 0) {
            Polygon best = null; double area = -1;
            for (int i=0;i<mp.getNumGeometries();i++) {
                Polygon p = (Polygon) mp.getGeometryN(i);
                if (p.getArea() > area) { area = p.getArea(); best = p; }
            }
            return best;
        }
        return sheet;
    }

    // aplica rotación y espejo a una geometría
    static Geometry orient(Geometry poly, int angleDeg, boolean mirrorX) {
        AffineTransformation at = new AffineTransformation();
        if (mirrorX) at = AffineTransformation.scaleInstance(-1, 1).compose(at);
        if (angleDeg % 360 != 0) {
            at = AffineTransformation.rotationInstance(Math.toRadians(angleDeg)).compose(at);
        }
        return at.transform(poly);
    }

    // traslada una geometría en dx y dy
    static Geometry translate(Geometry g, double dx, double dy) {
        return AffineTransformation.translationInstance(dx, dy).transform(g);
    }
}
