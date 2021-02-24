function tb = UnionBBoxes(tb)
tb.minpt = min(min(tb.bbox, [], 1), [], 3);
tb.maxpt = max(max(tb.bbox, [], 1), [], 3);