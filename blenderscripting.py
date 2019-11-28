#generate buildings with holes, the final script
import bpy
import numpy as np
import pyproj
from pykml import parser
from shapely.geometry import Polygon
from shapely.geometry import LinearRing
import shapely.geometry as sg


def split_horiz_by_point(polygon, point):
    """"""
    assert polygon.geom_type == "Polygon" and point.geom_type == "Point"
    nx, ny, xx, xy = polygon.bounds
    if point.x < nx or point.x > xx:
        return [polygon]
    print(point.x)
    lEnv = sg.LineString([(nx,      ny), (point.x, xy)]).envelope
    rEnv = sg.LineString([(point.x, ny), (xx,      xy)]).envelope
    
    try:
        return [polygon.intersection(lEnv), polygon.intersection(rEnv)]
    except Exception as e:
        print("Geometry error: %s" % validation.explain_validity(polygon))
        return [polygon.buffer(0)]


def check_split_multipoly(shape):
    """"""
    parts=[]
    if shape.type == "MultiPolygon":
        for p in shape.geoms:
            parts.extend(check_split_multipoly(p))
    elif shape.type == "Polygon":
        print(len(shape.interiors))
        if len(shape.interiors):
            pt = shape.interiors[0].centroid
            halves = split_horiz_by_point(shape, pt)
            print(halves)
            for p in halves:
                if p.type=="Polygon":
                    print('p',sg.mapping(p))
                    parts.extend(check_split_multipoly(p))
                    print("len of parts",len(parts))
                else:
                    for h in p:
                        if h.type=="Polygon":
                            print("h",sg.mapping(h))
                            parts.extend(check_split_multipoly(h))
        else:
            print('no interiors')
            parts = [shape]
    return parts




#filename='G:/data/thesisdataset/area10/KML/OGCKML_2/3880_5817_LoD2/doc.kml'
for i in range(13,17):
    if i<10:
        temp=str(0)+str(i)
    else:
        temp=str(i)
    name1='area'+temp
    print(name1)
    filename='/home/ge56cur/nas/ge56cur/masterthesis/KML/'+name1+'/'+name1+'.kml'
    cordiname='/home/ge56cur/nas/ge56cur/masterthesis/dataset1-linear/area'+temp+'/parameter/coordi.txt'
#filename='/home/ge56cur/nas/ge56cur/masterthesis/KML/area06/area06.kml' # to read kml
#cordiname='/home/ge56cur/nas/ge56cur/masterthesis/dataset1-linear/area06/parameter/coordi.txt'
    inProj = pyproj.Proj(init='epsg:4326') #WGS
    inProj1=pyproj.Proj(init='epsg:4978') #WGS geocentric   X,Y,Z
    outProj = pyproj.Proj(init='epsg:3857') #Web mercator
    outProj1 = pyproj.Proj(init='epsg:25833') #UTM 33N
    cordi=open(cordiname,"r")
    cordis=cordi.readlines()
    COORD=[]
    COORD.append(np.double(cordis[0].split(' ')))
    print(COORD)
    with open(filename) as f:
        folder = parser.parse(f).getroot().Document.Folder
        

    plnm=[]
    count=0
    cordi=[]
    cordi1=[]
    #for buildings in folder:
        #b_name=buildings.name
        #for wall in buildings.Placemark:
            #name=wall.name
            #print(name)
            #for poly in wall.MultiGeometry.Polygon:
                ##count=count+1
                #coord=poly.outerBoundaryIs.LinearRing.coordinates
                #plnm.append(name.text)
                #A=coord.text.split(' ')
                #cordi.append(A)
    for buildings in folder.Placemark:
        #name=buildings.name
        for poly in buildings.Polygon:
            count=count+1
            coord=poly.outerBoundaryIs.LinearRing.coordinates
            AA2=[]
            if hasattr(poly, 'innerBoundaryIs'):
                
                for L in poly.innerBoundaryIs:
                    
                    coord1=L.LinearRing.coordinates
                    AA=coord1.text.split(' ')
                    AA2.append(AA)
                print(len(AA2),count)
            #plnm.append(name.text)
            #bname.append(b_name.text)
            A=coord.text.split(' ')
            cordi.append(A)
            cordi1.append(AA2)
    print(len(cordi))

    ENU=[]
    for i in range(len(cordi)): #each polygon
        vetexx=[]
        vetexy=[]
        vetexz=[]
        for j in range(len(cordi[i])): #each vertex
            
            ss=cordi[i][j].split(',')
            #ENU fails
            #cc=geodetic_to_enu(float(ss[1]),float(ss[0]),float(ss[2]),camera_coord[0][0],camera_coord[0][1],camera_coord[0][2])
            #first longitude then latitude
            
            #x1,y1 = -11705274.6374,4826473.6922
            vetexx.append(float(ss[0]))
            vetexy.append(float(ss[1]))
            vetexz.append(float(ss[2]))
        vetexx=tuple(vetexx)
        vetexy=tuple(vetexy)
        vetexz=tuple(vetexz)
        #print(vetexx)
        cc=pyproj.transform(inProj,outProj1,vetexx,vetexy,vetexz)
        #print(cc[0]-390000,cc[1]-5819000)
        #zpp=zip(list(cc[0]*np.ones_like(cc[0])/100-3889.151999449386),list(cc[1]*np.ones_like(cc[1])/100-58179.26052129787),list(cc[2]*np.ones_like(cc[2])/100-69.66072412/100))
        zpp=zip(list(cc[0]*np.ones_like(cc[0])-COORD[0][0]),list(cc[1]*np.ones_like(cc[1])-COORD[0][1]),list(cc[2]*np.ones_like(cc[2])-COORD[0][2]-3.4))
        ENU.append(zpp)

    ENUin=[]
    for i in range(len(cordi1)):
        if len(cordi1[i]) == 0:
            ENUin.append(0)
        
        else:

            g=[]
            for j in range(len(cordi1[i])): #each group
                vetexx=[]
                vetexy=[]
                vetexz=[]
                for q in range(len(cordi1[i][j])): #each vetex
                    ss=cordi1[i][j][q].split(',')
                #ENU fails
                #cc=geodetic_to_enu(float(ss[1]),float(ss[0]),float(ss[2]),camera_coord[0][0],camera_coord[0][1],camera_coord[0][2])
                #first longitude then latitude
                
                #x1,y1 = -11705274.6374,4826473.6922
                    vetexx.append(float(ss[0]))
                    vetexy.append(float(ss[1]))
                    vetexz.append(float(ss[2]))
                cc=pyproj.transform(inProj,outProj1,vetexx,vetexy,vetexz)
                zpp=zip(list(cc[0]*np.ones_like(cc[0])-COORD[0][0]),list(cc[1]*np.ones_like(cc[1])-COORD[0][1]),list(cc[2]*np.ones_like(cc[2])-COORD[0][2]-3.4))
                g.append(zpp)
            ENUin.append(g)
            
            
            
     

    print(len(ENU))
    verts=[]
    faces=[]
    face=[]
    start=0

    for i in range(len(ENU)):
        p=ENU[i]
        p=list(p)
        ENU[i]=p
        p1=ENUin[i]
        inners=[]
        if p1!=0:
            p1_=[]
            for p11 in p1:
                p12=list(p11)
                p1_.append(p12)
                inners.append(Polygon(p12))
            ENUin[i]=p1_
        p1=ENUin[i]
        
        outer=Polygon(p)
        if p1!=0:
            print('len p',len(p),'len p1',len(p1))
            shape=Polygon(outer.exterior.coords,[inner.exterior.coords for inner in inners])
            parts1 = check_split_multipoly(shape)
            for pa in parts1:
                face=[]
                A=sg.mapping(pa)['coordinates']
                pb=list(A[0])
                verts=verts+pb
                face.append(tuple(list(range(start,start+len(pb)))))
                start=start+len(pb)
                faces=faces+face
        else:
            face=[]
            verts=verts+p
            face.append(tuple(list(range(start,start+len(p)))))
            start=start+len(p)
            faces=faces+face
            
    print(faces)

    print(faces)
    mymesh=bpy.data.meshes.new(name1)
    myobject=bpy.data.objects.new(name1,mymesh)
    #myobject.location=bpy.context.scene.cursor.location

    bpy.context.collection.objects.link(myobject)
    mymesh.from_pydata(verts,[],faces)
     
    mymesh.update(calc_edges=True)
        
    
    

