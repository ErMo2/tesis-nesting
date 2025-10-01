package org.tesis.nesting;

import org.locationtech.jts.geom.*;
import java.util.Locale;

class SvgWriter {

    static String toSVG(Polygon sheet, java.util.List<PlacedPart> placed) {
        Envelope env = sheet.getEnvelopeInternal();
        double minX = env.getMinX(), minY = env.getMinY(), w = env.getWidth(), h = env.getHeight();

        StringBuilder sb = new StringBuilder();
        sb.append("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        sb.append("<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"").append(w)
          .append("\" height=\"").append(h).append("\" viewBox=\"")
          .append(minX).append(" ").append(minY).append(" ").append(w).append(" ").append(h).append("\">\n");

        // Fondo
        sb.append("  <rect x=\"").append(minX).append("\" y=\"").append(minY)
          .append("\" width=\"").append(w).append("\" height=\"").append(h)
          .append("\" fill=\"white\"/>\n");

        // Grupo invertido en Y (SVG tiene Y hacia abajo)
        sb.append("  <g transform=\"translate(0,").append(minY + h).append(") scale(1,-1)\">\n");

        // Sheet
        sb.append("    <path d=\"").append(pathFor(sheet))
          .append("\" fill=\"#f0f0f0\" stroke=\"#333\" stroke-width=\"1\"/>\n");

        // Piezas
        String[] fills = {"#6baed6","#74c476","#fd8d3c","#9e9ac8","#fdd0a2","#a1d99b","#9ecae1","#fdae6b"};
        int colorIdx = 0;
        for (PlacedPart p : placed) {
            String fill = fills[colorIdx++ % fills.length];
            emitGeometryAsSvgPaths(sb, p.polyPlaced, fill, p);
        }

        sb.append("  </g>\n</svg>\n");
        return sb.toString();
    }

    // ---------- helpers de dibujo ----------

    static void emitGeometryAsSvgPaths(StringBuilder sb, Geometry g, String fill, PlacedPart meta) {
        if (g instanceof Polygon) {
            emitPolygon(sb, (Polygon) g, fill, meta);
        } else if (g instanceof MultiPolygon) {
            MultiPolygon mp = (MultiPolygon) g;
            for (int i = 0; i < mp.getNumGeometries(); i++) {
                emitPolygon(sb, (Polygon) mp.getGeometryN(i), fill, meta);
            }
        } else if (g instanceof GeometryCollection) {
            GeometryCollection gc = (GeometryCollection) g;
            for (int i = 0; i < gc.getNumGeometries(); i++) {
                emitGeometryAsSvgPaths(sb, gc.getGeometryN(i), fill, meta);
            }
        } else if (g instanceof LineString) {
            LineString ls = (LineString) g;
            sb.append("    <path d=\"");
            appendLineString(sb, ls);
            sb.append("\" fill=\"none\" stroke=\"#111\" stroke-width=\"0.8\"/>\n");
        }
    }

    static void emitPolygon(StringBuilder sb, Polygon poly, String fill, PlacedPart meta) {
        sb.append("    <path d=\"").append(pathFor(poly))
          .append("\" fill=\"").append(fill)
          .append("\" fill-opacity=\"0.85\" stroke=\"#111\" stroke-width=\"0.8\">\n");
        if (meta != null) {
            sb.append("      <title>").append(meta.id)
              .append(" | θ=").append(meta.angleDeg)
              .append("° | mirror=").append(meta.mirrored)
              .append("</title>\n");
        }
        sb.append("    </path>\n");
    }

    // Genera el atributo "d" de un path SVG a partir de un Polygon (exterior + huecos)
    static String pathFor(Polygon poly) {
        StringBuilder sb = new StringBuilder();
        appendLineString(sb, poly.getExteriorRing());
        for (int i = 0; i < poly.getNumInteriorRing(); i++) {
            appendLineString(sb, poly.getInteriorRingN(i));
        }
        return sb.toString();
    }

    // Agrega comandos M/L/Z para una LineString cerrada
    static void appendLineString(StringBuilder sb, LineString ls) {
        Coordinate[] c = ls.getCoordinates();
        if (c.length == 0) return;
        sb.append("M ").append(fmt(c[0].x)).append(" ").append(fmt(c[0].y)).append(" ");
        for (int i = 1; i < c.length; i++) {
            sb.append("L ").append(fmt(c[i].x)).append(" ").append(fmt(c[i].y)).append(" ");
        }
        sb.append("Z ");
    }

    static String fmt(double d) {
        return String.format(Locale.US, "%.3f", d);
    }
}
