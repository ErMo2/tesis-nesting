package org.tesis.nesting;

import org.locationtech.jts.geom.Geometry;

public class PlacedPart {
    String   id;
    Geometry polyPlaced;   // Polygon / MultiPolygon (por buffers/transformaciones)
    double   angleDeg;
    boolean  mirrored;
    double   dx, dy;
}
