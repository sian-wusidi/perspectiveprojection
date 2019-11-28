# perspectiveprojection
to project a 3D city model (KML) under along a camera trajectory
1. parameter_project_final.py is the original implementation in python, the occlusion is approximated to render with the following sequence: walls->roofs (outer polygons)-> roofs (inner polygons)
2. blenderscripting.py is the script in blender to read and render 3D city model from KML, including the split of rings into small polygons
