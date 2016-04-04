def get_line_intersection( p0_x,  p0_y,  p1_x,  p1_y, p2_x,  p2_y,  p3_x,  p3_y):
    s10_x = p1_x - p0_x
    s10_y = p1_y - p0_y
    s32_x = p3_x - p2_x
    s32_y = p3_y - p2_y

    denom = s10_x * s32_y - s32_x * s10_y
    if denom == 0:
        return False # Collinear
    denomPositive = denom > 0

    s02_x = p0_x - p2_x
    s02_y = p0_y - p2_y
    s_numer = s10_x * s02_y - s10_y * s02_x
    if (s_numer < 0) == denomPositive:
        return False # No collision

    t_numer = s32_x * s02_y - s32_y * s02_x
    if (t_numer < 0) == denomPositive:
        return False # No collision

    if ((s_numer > denom) == denomPositive) or ((t_numer > denom) == denomPositive):
        return False 
    # Collision detected
    t = t_numer / denom
    i_x = p0_x + (t * s10_x)
    i_y = p0_y + (t * s10_y)

    return (True,(i_x,i_y))

print get_line_intersection(0,0,2,0, 1,1,1,-1)