package org.tesis.nesting;

import java.util.Locale;

public class VertexRow {
    String entityType;   // SHEET | PART
    String entityId;
    int ring;
    int pointIdx;
    double x, y;

    Integer qty;                 // PART
    Integer rotationStepDeg;     // PART
    Boolean mirrorAllowed;       // PART
    Double clearance;            // PART
    Double grainAngleDeg;        // PART (no usado aquí)

    static VertexRow fromCsv(String[] h, String[] v) {
        VertexRow r = new VertexRow();
        r.entityType = get(v, idx(h, "entity_type"));
        r.entityId   = get(v, idx(h, "entity_id"));
        r.ring       = parseInt(get(v, idx(h, "ring")), 0);
        r.pointIdx   = parseInt(get(v, idx(h, "point_idx")), 0);
        r.x          = parseDouble(get(v, idx(h, "x")), 0);
        r.y          = parseDouble(get(v, idx(h, "y")), 0);

        r.qty              = parseNullableInt(getOpt(v, h, "qty"));
        r.rotationStepDeg  = parseNullableInt(getOpt(v, h, "rotation_step_deg"));
        r.mirrorAllowed    = parseNullableBool(getOpt(v, h, "mirror_allowed"));
        r.clearance        = parseNullableDouble(getOpt(v, h, "clearance"));
        r.grainAngleDeg    = parseNullableDouble(getOpt(v, h, "grain_angle_deg"));
        return r;
    }

    // ------------- helpers CSV -------------
    static String[] splitCsv(String s) {
        String[] raw = s.split(",", -1);
        for (int i=0;i<raw.length;i++) raw[i] = raw[i].trim();
        return raw;
    }
    static int idx(String[] h, String name) {
        for (int i=0;i<h.length;i++) if (h[i].equalsIgnoreCase(name)) return i;
        throw new IllegalArgumentException("Cabecera CSV faltante: " + name);
    }
    static int idxOpt(String[] h, String name) {
        for (int i=0;i<h.length;i++) if (h[i].equalsIgnoreCase(name)) return i;
        return -1;
    }
    static String get(String[] v, int idx) {
        if (idx < 0 || idx >= v.length) return "";
        return v[idx];
    }
    static String getOpt(String[] v, String[] h, String name) {
        int i = idxOpt(h, name);
        return i == -1 ? "" : get(v, i);
    }

    static Integer parseNullableInt(String s) {
        if (s == null || s.isEmpty()) return null;
        try { return Integer.parseInt(s); } catch (Exception e) { return null; }
    }
    static Double parseNullableDouble(String s) {
        if (s == null || s.isEmpty()) return null;
        try { return Double.parseDouble(s); } catch (Exception e) { return null; }
    }
    static Boolean parseNullableBool(String s) {
        if (s == null || s.isEmpty()) return null;
        String t = s.toLowerCase(Locale.ROOT);
        if (t.equals("true") || t.equals("1") || t.equals("yes") || t.equals("y")) return true;
        if (t.equals("false") || t.equals("0") || t.equals("no") || t.equals("n")) return false;
        return null;
    }
    static int parseInt(String s, int d) {
        Integer v = parseNullableInt(s); return v == null ? d : v;
    }
    static double parseDouble(String s, double d) {
        Double v = parseNullableDouble(s); return v == null ? d : v;
    }
}
