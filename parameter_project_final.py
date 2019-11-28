# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 14:00:25 2019

@author: Sidi Wu
"""
from pykml import parser
import pandas as pd
import pyproj
import numpy as np
import math
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import json
from coordinate_transfer import ecef_to_enu,ECEF_to_geodetic
import cv2
import os
import copy



def takeThird(elem):
    return elem[2]  
def calcangle(a1,a2):
    
    if len(np.shape(a1))>1 or len(np.shape(a2))>1:
        #print('true')
        if np.shape(a1)[0]==3:
            a1=np.transpose(a1)
        if np.shape(a2)[1]==3:
            a2=np.transpose(a2)
    
    b=np.dot(a1,a2)/(np.linalg.norm(a1)*np.linalg.norm(a2))
    #print(b>=1.0)
    if b>=1.0:
        a=math.degrees(math.acos(1))
    elif b<=-1.0:
        a=math.degrees(math.acos(-1))
    else:
        a=math.degrees(math.acos(b))
    
    return a
def calcrotation(alpha,beta,gamma):
    alpha=math.radians(alpha)
    beta=math.radians(beta)
    gamma=math.radians(gamma)
    #print(gamma)
    R=np.zeros([3,3])
    R1=np.mat([[1,0,0],[0,math.cos(alpha),math.sin(alpha)],[0,-math.sin(alpha),math.cos(alpha)]])
    R2=np.mat([[math.cos(beta),0,-math.sin(beta)],[0,1,0],[math.sin(beta),0,math.cos(beta)]])
    R3=np.mat([[math.cos(gamma),math.sin(gamma),0],[-math.sin(gamma),math.cos(gamma),0],[0,0,1]])
    #print(R3)
    #y(roll)-x(pitch/tilt)-z(yaw/pan) the multiplication sequence doenst matter much
    #R=np.matmul(R2,np.matmul(R1,R3))
    R=np.matmul(R1,np.matmul(R2,R3))
    #print(R)
    #print(R3)
    return R

def projection2d(R,T,K,X):
    Xc=np.ones([4,1])
    Xc[0:3,:]=X+T
    
    #print(X+T)
    #point to object
    if (Xc[2,:]>0):
        print(Xc[2,:])
    Xc[2,:]=np.abs(Xc[2,:])
    Xc[0:3,:]=np.dot(R,Xc[0:3,:])
    Xc[3]=1
    #translation already contained in the calibration matrix!!!
    Xp=np.dot(K,Xc)
    Xp=Xp/Xp[2]
    return Xp,Xc

        
        
def equation_plane(p1,p2,p3):  
    x1=p1[0]
    y1=p1[1]
    z1=p1[2]
    x2=p2[0]
    y2=p2[1]
    z2=p2[2]
    x3=p3[0]
    y3=p3[1]
    z3=p3[2]
    a1 = x2 - x1 
    b1 = y2 - y1 
    c1 = z2 - z1 
    a2 = x3 - x1 
    b2 = y3 - y1 
    c2 = z3 - z1 
    a = b1 * c2 - b2 * c1 
    b = a2 * c1 - a1 * c2 
    c = a1 * b2 - b1 * a2 
    d = (- a * x1 - b * y1 - c * z1) 
    #print( "equation of plane is "), 
    #print (a, "x +"), 
    #print (b, "y +"), 
    #print (c, "z +"), 
    #print (d, "= 0.")
    
    return(a,b,c,d)
    
def calcequaplane(P):
    SUM=[]
    PLANE=[]
    for i in range(0,len(P)-2):
        sum1=0
        plane=equation_plane(P[i],P[i+1],P[i+2])
        PLANE.append(plane)
    #index=random.sample(range(0,len(P)), 3)
#    plane=equation_plane(P[index[0]],P[index[1]],P[index[2]])
        for p in P:
            sum1=sum1+plane[0]*p[0]+plane[1]*p[1]+plane[2]*p[2]+plane[3]
        SUM.append(sum1)
    
    pos=SUM.index(min(SUM))
    return PLANE[pos],min(SUM)

inProj = pyproj.Proj(init='epsg:4326') #WGS
inProj1=pyproj.Proj(init='epsg:4978') #WGS geocentric   X,Y,Z
outProj = pyproj.Proj(init='epsg:3857') #Web mercator
outProj1 = pyproj.Proj(init='epsg:25833') #UTM 33N

#%% read file and get the pan and roll
dics1=os.listdir()
for dic in dics1:
    if os.path.isdir(dic)==True:
        print(dic)
dic_=dics1[10]
if os.path.isdir(dic_)==True: 
    os.chdir(dic_)
CORC=np.array([[1,0,0],[0,1,0],[0,0,1]])
PAN=[]
PITCH=[]
XYZ=[]
Height=[]
CROSS=[]
file1 = open("metadata.txt","r") 
lines = file1.read().split(',')
groundheight=float(lines[1])
heightstart=float(lines[3])
Height_wrt_ground=[]
with open(dic_+'.json') as json_file:  # 需要等于文件名
    data = json.load(json_file)
for p in data['cameraFrames']:
    xyz=np.array([p['position']['x'],p['position']['y'],p['position']['z']])
    XYZ.append(xyz)
    R=calcrotation(180-p['rotation']['x'],p['rotation']['y'],p['rotation']['z'])
    TEMP=np.dot(np.transpose(R),CORC)

    E_x=ecef_to_enu(TEMP[0,0]+xyz[0],TEMP[0,1]+xyz[1],TEMP[0,2]+xyz[2], xyz[0], xyz[1], xyz[2])
    N_y=ecef_to_enu(TEMP[1,0]+xyz[0],TEMP[1,1]+xyz[1],TEMP[1,2]+xyz[2], xyz[0], xyz[1], xyz[2])
    U_z=ecef_to_enu(TEMP[2,0]+xyz[0],TEMP[2,1]+xyz[1],TEMP[2,2]+xyz[2], xyz[0], xyz[1], xyz[2])
# 算出来y,z都是旋转了40°, y的旋转有pan，也有tilt，所以只考虑x的夹角-pan, z的夹角-roll
    #CROSS.append(np.cross(np.array([[1,0,0]]),np.array([E_x])))
    pan=calcangle(np.array([[1,0,0]]),np.array([E_x]))
    if np.cross(np.array([[1,0,0]]),np.array([E_x]))[0,2]>0:
        pan=-pan
    pitch=calcangle(np.array([[0,0,1]]),np.array([U_z]))
    PAN.append(pan)
    PITCH.append(pitch)
    camera_coord=ECEF_to_geodetic(xyz[0],xyz[1],xyz[2])
    camera_proj=[pyproj.transform(inProj,outProj1,camera_coord[1],camera_coord[0],camera_coord[2])]
    
    Height.append(camera_proj[0][2])
    
for i in range(len(Height)):
    Height_wrt_ground.append(heightstart+Height[i]-Height[0])

os.mkdir('parameter')  
np.savetxt('parameter/pan.txt', PAN, delimiter=' ')   # X is an array
np.savetxt('parameter/pitch.txt', PITCH, delimiter=' ')   # X is an array
np.savetxt('parameter/height.txt', Height_wrt_ground, delimiter=' ')   # X is an array

#K=np.array([[1153.88428,0,959,0],[0,1153.88428,539,0],[0,0,1,0]])  #optimization for ETRS fov=50
K=np.array([[768.4365663420178,0,959,0],[0,768.4365663420178,539,0],[0,0,1,0]]) #fov=70
f2=768.4365663420178
R=calcrotation(-PITCH[35],0,-PAN[35]-1.1) #fov=50 -1.5
T=np.zeros([3,1])

#%% READ KML
path='KML/OGCKML_2'
dics=os.listdir(path)
count=0
plnm=[]
cordi=[]
cordi1=[]
bname=[]
for dic in dics:
    if os.path.isdir(path+'/'+dic)==True:
        filename=path+'/'+dic+'/doc.kml'
        count=count+1
        print(filename)
#filename='doc.kml'
        with open(filename) as f:
            folder = parser.parse(f).getroot().Document.Folder
    
        for buildings in folder:
            AA=[]
            
            b_name=buildings.name
            for wall in buildings.Placemark:
                name=wall.name
                for poly in wall.MultiGeometry.Polygon:
                    count=count+1
                    coord=poly.outerBoundaryIs.LinearRing.coordinates
                    AA2=[]
                    if hasattr(poly, 'innerBoundaryIs'):
                        
                        for L in poly.innerBoundaryIs:
                            
                            coord1=L.LinearRing.coordinates
                            AA=coord1.text.split(' ')
                            AA2.append(AA)
                        print(len(AA2),count)
                    plnm.append(name.text)
                    bname.append(b_name.text)
                    A=coord.text.split(' ')
                    cordi.append(A)
                    cordi1.append(AA2)




# accelerate by zip first!!
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
    cc=pyproj.transform(inProj,outProj1,vetexx,vetexy,vetexz)
    zpp=zip(list(cc[0]),list(cc[1]),list(cc[2]))
    ENU.append(zpp)
# accelerate by zip first!!
ENUin=[]
for i in range(len(ENU)):
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
            zpp=zip(list(cc[0]),list(cc[1]),list(cc[2]))
            g.append(zpp)
        ENUin.append(g)


#%% prepare ground truth
# folder for semantics
os.mkdir('semantics')  
os.mkdir('depth')
os.mkdir('edges')  
os.mkdir('mask')  
os.mkdir('surfaces')    
os.mkdir('semanticscheck')    
point_size = 20
point_color = (0, 0, 255, 200) # BGR
point_color1 = (0, 255, 0, 200) # BGR
thickness = 4 # 可以为 0 、4、8
#project2D=np.zeros([1080,1920])




for idx in range(len(XYZ)):
    camera_coord=ECEF_to_geodetic(XYZ[idx][0],XYZ[idx][1],XYZ[idx][2])
    camera_proj=[pyproj.transform(inProj,outProj1,camera_coord[1],camera_coord[0],camera_coord[2])]
    R=calcrotation(-PITCH[idx],0,-PAN[idx]-1.5)
    T[:,0]=np.array([-camera_proj[0][0],-camera_proj[0][1],-Height_wrt_ground[idx]+groundheight]) #after optimization
    count=0
    COL=[]
    ROW=[]
    temp=[]
    #open cv fill the polygons
    projectpoly1=[] #roof
    #projectpoly1in=[] #roof inner ring
    projectpoly2=[] #wall
    #projectpoly2in=[] #wall inner ring
    # please pay attention to the height!!!
   
    for i in range(len(ENU)):
        sum1=0
        p=ENU[i]
        p1=ENUin[i]
        p=list(p)
        if p1!=0:
            p1_=[]
            for p11 in p1:
                p12=list(p11)
                p1_.append(p12)
            ENUin[count]=p1_
        p1=ENUin[i]
        ENU[count]=p
        point_list1=[]
        point_list2=[]
        point_list1ccs=[]
        point_list2ccs=[]
        A=[]
        depth0=[]

        if plnm[count]=='RoofSurface':
            
            
            for i in range(0,len(p)):
                B=np.zeros([3,1])
                B[:,0]=np.transpose(p[i])
                B1,B2=projection2d(R,T,K,B)
                col=math.floor(B1[0,0])
                row=math.floor(B1[1,0])
                point_list1.append([col,1079-row])
                A.append(B2)
                #depth0.append(B2[2][0])
                depth0.append(np.sqrt(pow(B2[2][0],2)+pow(B2[1][0],2)+pow(B2[0][0],2)))
                #contours.append((col,1079-row))
                #print((B+T),'and',B1)
            d=np.max(depth0)
            poly1=np.asarray(point_list1)
            #LA=len(A)
            #plane=equation_plane(A[np.int((LA+1)/2)-1],A[np.int((LA+1)/2)],A[np.int((LA+1)/2)+1])
            #con=np.asarray(contours)
            
            if p1!=0:
                
                for i in range(0,len(p1)):
                    point_list1ccs1=[]
                    for j in range(0,len(p1[i])):
                        p12=list(p1[i])
                        B=np.zeros([3,1])
                        B[:,0]=np.transpose(p12[j])
                        B1,B2=projection2d(R,T,K,B)
                        col=math.floor(B1[0,0])
                        row=math.floor(B1[1,0])
                        point_list1ccs1.append([col,1079-row])
                    point_list1ccs.append(point_list1ccs1)
                print(count)
                
                print(len(point_list1ccs))
            poly11=np.asarray(np.asarray(point_list1ccs))
            
            
            if (poly1[:,0]<1920).any() and (poly1[:,1]<1079).any() and (poly1[:,0]>=0).any() and (poly1[:,1]>=0).any():
    
                plane,sum1=calcequaplane(A)
                projectpoly1.append([poly1,plane,d,poly11]) # the last is the inside
                #projectpoly1in.append(poly11)
            #contourpoly.append(con)
        #if (data['features'][count]['attributes']['surface']=='wall'):
        if plnm[count]=='WallSurface':
            for i in range(0,len(p)):
                B=np.zeros([3,1])
                B[:,0]=np.transpose(p[i])
                B1,B2=projection2d(R,T,K,B)
                col=math.floor(B1[0,0])
                row=math.floor(B1[1,0])
                point_list2.append([col,1079-row])
                A.append(B2)
                #depth0.append(B2[2][0])
                depth0.append(np.sqrt(pow(B2[2][0],2)+pow(B2[1][0],2)+pow(B2[0][0],2)))
                #contours.append((col,1079-row))
                #print((B+T),'and',B1)
            d=np.max(depth0)
            poly2=np.asarray(point_list2)
            #LA=len(A)
            #plane=equation_plane(A[np.int((LA+1)/2)-1],A[np.int((LA+1)/2)],A[np.int((LA+1)/2)+1])
            if p1!=0:
                
                point_list1ccs1=[]
                for i in range(0,len(p1)):
                    for j in range(0,len(p1[i])):
                        p12=list(p1[i])
                        B=np.zeros([3,1])
                        B[:,0]=np.transpose(p12[j])
                        B1,B2=projection2d(R,T,K,B)
                        col=math.floor(B1[0,0])
                        row=math.floor(B1[1,0])
                        point_list1ccs1.append([col,1079-row])
                point_list2ccs.append(point_list1ccs1)
            poly22=np.asarray(np.asarray(point_list2ccs))
            

            
            
            if (poly2[:,0]<1920).any() and (poly2[:,1]<1079).any() and (poly2[:,0]>=0).any() and (poly2[:,1]>=0).any():
    
                plane,sum1=calcequaplane(A)
                projectpoly2.append([poly2,plane,d,poly22])
                #projectpoly2in.append(poly22)
                
                
                
            if sum1>10:
                print(sum1)
            #con=np.asarray(contours)
            
            #projectpolyplane2.append(plane)
    #        #contourpoly.append(con)
        count=count+1
        
    projectpoly1.sort(key=takeThird,reverse = True)
    projectpoly2.sort(key=takeThird,reverse = True)
   # create labels
    img = np.zeros((1080, 1920, 3), np.uint8)
    img2= np.zeros((1080, 1920, 3), np.uint8)
    A2=np.zeros((1080, 1920), np.uint8)
    depth=np.zeros((1080, 1920, 1),np.double)
    label=np.zeros((1080, 1920, 3),np.uint8)
    surface=np.zeros((1080, 1920, 3),np.uint8)
    for i in range(0,len(projectpoly2)):
        tt=[] 
        img = np.zeros((1080, 1920, 3), np.uint8)
        polygon=projectpoly2[i][0]
        plane=projectpoly2[i][1]
        polygonin=projectpoly2[i][3]
        increpoly=[]
        color1=np.random.randint(50,255)
    
        tt=cv2.fillConvexPoly(img,polygon,point_color1)
        inn=np.ones_like(depth)
        inn1=np.ones_like(label)       
        increpoly=np.where((tt[:,:,1]) !=0)
    
        for j in range(0,np.shape(increpoly)[1]):#调用的时候先行后列
            x=increpoly[1][j]-959
            y=1079-increpoly[0][j]-539
            
            pl=-plane[3]*f2/(plane[0]*x+plane[1]*y+plane[2]*f2)
#            X=a*x/f2
#            Y=a*y/f2
#            d=np.sqrt(pow(X,2)+pow(Y,2)+pow(a,2))
#            X1=depth[increpoly[0][j],increpoly[1][j]]*x/f2
#            Y1=depth[increpoly[0][j],increpoly[1][j]]*y/f2
#            d1=np.sqrt(pow(X1,2)+pow(Y1,2)+pow(depth[increpoly[0][j],increpoly[1][j]],2))

            
            depth[increpoly[0][j],increpoly[1][j]]=pl
            label[increpoly[0][j],increpoly[1][j]]=(0,0,255)
            surface[increpoly[0][j],increpoly[1][j]]=(0,0,color1)
            
        for q in polygonin:
            img = np.zeros((1080, 1920, 3), np.uint8)

            tt=[]
            tt=cv2.fillConvexPoly(img,np.asarray(q),point_color1)
            increpoly=np.where((tt[:,:,1]) !=0)
            inn[increpoly]=0
            inn1[increpoly]=(0,0,0)
    
    depth=depth*inn
    label=label*inn1
    surface=surface*inn1
            


    

        
    depth11=copy.copy(depth)
    for i in range(0,len(projectpoly1)):
        polygon=projectpoly1[i][0]
        plane=projectpoly1[i][1]
        polygonin=projectpoly1[i][3]

        
        color1=np.random.randint(50,255)
  
        inn=np.ones_like(depth)
        inn1=np.ones_like(label)
        for q in polygonin:
            img = np.zeros((1080, 1920, 3), np.uint8)

            tt=[]
            tt=cv2.fillConvexPoly(img,np.asarray(q),point_color1)
            increpoly=np.where((tt[:,:,1]) !=0)
            inn[increpoly]=0
            inn1[increpoly]=(0,0,0)
        
        tt=[] 
        tt=cv2.fillConvexPoly(img,polygon,point_color1)
        increpoly=[]
        increpoly=np.where((tt[:,:,1]) !=0)
        img = np.zeros((1080, 1920, 3), np.uint8)

        for j in range(0,np.shape(increpoly)[1]):#调用的时候先行后列
            x=increpoly[1][j]-959
            y=1079-increpoly[0][j]-539
            
            pl=-plane[3]*f2/(plane[0]*x+plane[1]*y+plane[2]*f2)
#            X=a*x/f2
#            Y=a*y/f2
#            d=np.sqrt(pow(X,2)+pow(Y,2)+pow(a,2))
#            X1=depth[increpoly[0][j],increpoly[1][j]]*x/f2
#            Y1=depth[increpoly[0][j],increpoly[1][j]]*y/f2
#            d1=np.sqrt(pow(X1,2)+pow(Y1,2)+pow(depth[increpoly[0][j],increpoly[1][j]],2))

            
            if pl*inn[increpoly[0][j],increpoly[1][j]]>0 and ((depth[increpoly[0][j],increpoly[1][j]]-pl*inn[increpoly[0][j],increpoly[1][j]]>=-0.2) or (depth[increpoly[0][j],increpoly[1][j]]==0)):
                depth11[increpoly[0][j],increpoly[1][j]]=pl*inn[increpoly[0][j],increpoly[1][j]]
                label[increpoly[0][j],increpoly[1][j]]=(0,0,255)*inn1[increpoly[0][j],increpoly[1][j]]
                surface[increpoly[0][j],increpoly[1][j]]=(0,0,color1) *inn1[increpoly[0][j],increpoly[1][j]]
                
                
    depth1=depth11
    temp=np.nonzero(depth1[:,:,0])
    depth1[temp]=(depth1[temp]-np.min(depth1[temp]))/(np.max(depth1[temp])-np.min(depth1[temp]))*255
    depth1 = cv2.resize(np.uint8(depth1), (480, 270), interpolation=cv2.INTER_CUBIC)
    label= cv2.resize(label, (480, 270), interpolation=cv2.INTER_CUBIC)
    surface= cv2.resize(surface, (480, 270), interpolation=cv2.INTER_CUBIC)
    nn=np.where(label[:,:]!=np.array([0,0,0]))
    mask=np.zeros_like(label)
    mask[nn[0],nn[1]]=np.array([1,1,1])*255
    mask1=cv2.resize(mask,(480, 270), interpolation=cv2.INTER_CUBIC)   
    edges= cv2.Canny(surface,50,100)
         
    if idx<10:
        img1name='footage/'+dic_+'_00'+str(idx)+'.jpeg'
    elif idx<100:   
        img1name='footage/'+dic_+'_0'+str(idx)+'.jpeg'
    else:   
        img1name='footage/'+dic_+'_'+str(idx)+'.jpeg'    
    img1=cv2.imread(img1name)
#    img = np.zeros((1080, 1920, 3), np.uint8)
    #480,270
#    for polygon in projectpoly2:
#        cv2.fillConvexPoly(img,polygon,point_color1)
#    
#    for polygon in projectpoly1:
#        cv2.fillConvexPoly(img,polygon,point_color)
#    img = cv2.resize(img, (480,270), interpolation=cv2.INTER_CUBIC)
    alpha = 0.6
    beta = 1-alpha
    gamma = 0
    img_add = cv2.addWeighted(img1, alpha,label, beta, gamma)
    #cv2.namedWindow("new18")
    #cv2.imshow('new18', img_add)
    imgwrname='semanticscheck/'+str(idx)+'.jpg'
    print(imgwrname)
    cv2.imwrite(imgwrname,img_add)
    imgwrname='semantics/'+str(idx)+'.jpg'
    cv2.imwrite(imgwrname,label)
    imgwrname='depth/'+str(idx)+'.jpg'
    cv2.imwrite(imgwrname,depth1)
    imgwrname='mask/'+str(idx)+'.jpg'
    cv2.imwrite(imgwrname,mask1)
    imgwrname='surfaces/'+str(idx)+'.jpg'
    cv2.imwrite(imgwrname,surface)
    imgwrname='edges/'+str(idx)+'.jpg'
    cv2.imwrite(imgwrname,edges)    
        

    


