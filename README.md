# Visibility-Algorithm
## Demonstrations
### Light Fading Visualization
![](docs/VisibilityAlgorithmPretty-6_17_2022.mp4)
### Ray Intersection Visualization
![](https://github.com/Saccharine-Coal/Visibility-Algorithm/blob/72bed7021f91510d80b5972b005a96a640eafbea/docs/VisibilityAlgorithmPretty-6_17_2022.mp4)
## Description
This program utilizes segment segment intersections to calculate polygon-ray intersections from a given origin. With these intersections a simple  triangulation algorithm constructs triangles, or polygons, that represent the visibility relative to the oring. This visibilit is visualized with the Python library Pygame. Numpy is largely used for vector operations.
## Performance
Performance, at best, decreases linearly with the number of vertices present. If we measure performance with fps (frames per second) and let n be the number of vertices, then the performance of the algorithm can be defined as:
$$ performance(n) \propto {1 \over n} $$
In reality, the performance is much worse then this linear relationship. For example: 
- 5 vertices results in a steady 30 fps.
- 35 vertices results in 5-15 fps.
## Dependencies
The dependencies for this project are:
- Numpy
- Pygame 2
- Shapely
- Python math module
- Python 3.5 or greater (type annotations)
- Python 3.7 or greater (dataclasses)
## Further Reading
- ...
## TODO List
### Initial Commit
- [X] Upload working tree
- [X] Complete README outline
- [X] Add gifs of application in action
- [ ] Log current bugs in README
- [X] Add dependencies
### Initial Refactor
- [ ] Remove commented code
- [ ] Remove unused code
- [ ] Remove unused modules
- [ ] Move similar code in same modules
- [ ] Move similar modules in same libs
- [ ] Add function documentation
### Current Bug Squashing
- [ ] ...
### Future Plans
- [ ] Implemenent add, sub, mult, div dunder methods in Point class
- [ ] Develop robust concave point in polygon function to eliminate shapely dependency
