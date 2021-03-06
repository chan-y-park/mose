### Temporary knobs to try different combinations
### of the available algorithms

# Add branch points at the ends of a MS 
# wall generated by two primary K-walls
ENHANCE_MS_WALLS = True

# Remove odd points from an MS wall
# best if used with MS_WALL_SORTING
# set to 'sweep' or 'phase'
CLEANUP_MS_WALLS = True

# If false, will use charge orbits
SORT_BY_GENEALOGY = True

# If false, will use straight approximation
CUT_K_WALLS = True

# Three options: 'sweep', 'neighbor', 'phase'
MS_WALLS_SORTING = 'phase'

# Whether to connect intersection points of an MS wall with a line
PLOT_MS_WALL_LINKS = True

# Decide whether to ignore intersections with pairing >2 or <0
# will make the code run much faster
IGNORE_WILD_INTERSECTIONS = True

