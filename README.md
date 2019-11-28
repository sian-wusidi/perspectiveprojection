# perspectiveprojection
to project a 3D city model (KML) under along a camera trajectory
parameter_project_final.py is the original implementation in python, the occlusion is approximated to render with the following sequence: walls->roofs (outer polygons)-> roofs (inner polygons)
blenderscripting.py is the script in blender to read and render 3D city model from cityGML, including the split of rings into small polygons
