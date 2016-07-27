library(plotly)
library(htmlwidgets)
library(proxy)

x.bga<-bga(y0$E,classvec = TS)
sampleCoords<-x$bet$ls
#rownames(sampleCoords)<-TS
sampleCentroid<-by(sampleCoords, TS, function(x) c((1/3)*sum(x[,1]),(1/3)*sum(x[,2]),(1/3)*sum(x[,3])), simplify = F )
sampleCentroid<-do.call(rbind, sampleCentroid)
sampleRanges<-by(sampleCoords, TS, function(x) apply(x, 2, range))
colnames(sampleCentroid)<-c("CS1","CS2","CS3")
dat<-rbind(sampleCoords,sampleCentroid)
groups<-as.factor(c( as.character(TS),
                     rep("Centroid",4)))

s1<-sampleCoords[TS==TS[3],]
s1.dist<-dist(s1, t(sampleCentroid[TS[3],]), method="euclidean")
s1.range<-apply(s1, 2, function(x) diff(range(x)))
s1.range<-s1.range/max(s1.range)
theta<-seq(from=0,to=2*pi,length.out = 100)
phi<-seq(0, pi, length.out = 100)
x<-max(s1.dist)*cos(theta) %o% sin(phi)/s1.range[1]
y<-max(s1.dist)*sin(theta) %o% sin(phi)/s1.range[2]
z<-rep(max(s1.dist),100) %o% cos(phi)/s1.range[3]



theta1<-seq(from=0,to=2*pi,length.out = 100)
phi1<-seq(0, pi, length.out = 100)

x1<-cos(theta1) %o% sin(phi1)
y1<-sin(theta) %o% sin(phi)
z1<-rep(1, 100) %o% cos(phi)
x2<-x/200
y2<-y/200
z2<-z/200

p<-plot_ly(dat, x=CS1,y=CS2, z=CS3, 
        hoverinfo="text", type="scatter3d", 
        mode="markers",group = groups, text=groups, marker=list(size=c(rep(12,12), rep(6,4))))
p<-add_trace(p, type="surface", x=x+sampleCentroid[TS[3],1], y=y+sampleCentroid[TS[3],2], z=z+sampleCentroid[TS[3],3],
          opacity=0.98, showscale=FALSE)
# p<-add_trace(p, type="surface", x=x1+sampleCentroid[2,1], y=y1+sampleCentroid[2,2], z=z1+sampleCentroid[2,3],
#              opacity=0.98, showscale=FALSE)
# p<-add_trace(p, type="surface", x=x1+sampleCentroid[3,1], y=y1+sampleCentroid[3,2], z=z1+sampleCentroid[3,3],
#              opacity=0.98, showlegend=FALSE, showscale=FALSE)
# p<-add_trace(p, type="surface", x=x/100+sampleCentroid[4,1], y=y/35+sampleCentroid[4,2], 
#              z=z/150+sampleCentroid[4,3],
#              opacity=0.98, showlegend=FALSE, showscale=FALSE)
p

saveWidget(as.widget(p), "~/Desktop/bga_sample.html")



p1<-plot_ly(x=x1,y=y1,z=z1, type="surface",opacity=0.5)
add_trace(p1, x=x, y=y, z=z, type="surface")








fac.bga<-x$bet$ls
colnames(fac.bga)<-c("Comp1","Comp2","Comp3")
probe.bga<-x$ord$ord$co[,1:3]

probe.bga<-rbind(fac.bga,probe.bga)
groups<-as.factor(c( as.character(TS),
                     rep("Gene",length(y0$genes$SystematicName))))

plot_ly(data = probe.bga, x=Comp1,y=Comp2, z=Comp3, group = groups,
        type="scatter3d", 
        mode="markers",
        text=c(as.character(TS),y0$genes$SystematicName),
        #marker=list(size=c(rep(12,12),rep(6,41000))),
        hoverinfo="text+group"
        )


topgenes(x.bga,axis = 1,n = 10)
topgenes(x.bga,axis = 2,n = 10)
topgenes(x.bga,axis = 3,n = 10)


# There are a number of ways of characterizing a three-dimensional rotation
# about the origin. The way I find easiest to remember is this. Let r = [u;v;w]
# be a unit vector which points along the axis of desired rotation and let rho be
# the desired angle of right-hand rotation to be made about this axis. Let p =
#     [x;y;z] be the coordinates of an arbitrary point. Then the corresonding
# rotated point P = [X;Y;Z] will be:
#     
#     P = dot(r,p)*r + cos(rho)*cross(cross(r,p),r) + sin(rho)*cross(r,p);
# 
# This can be rewritten in terms of a fixed 3 x 3 rotation matrix, R, to be
# multiplied by any p:
#     
#     cr = cos(rho); sr = sin(rho);
# R = [(1-cr)*u^2+cr,(1-cr)*u*v-sr*w,(1-cr)*w*u+sr*v;
#      (1-cr)*u*v+sr*w,(1-cr)*v^2+cr,(1-cr)*v*w-sr*u;
#      (1-cr)*w*u-sr*v,(1-cr)*v*w+sr*u,(1-cr)*w^2+cr];
# 
# P = R*p;
# 
# Then you can make the desired translation. Of course in your case you will
# apply this to the points of the parameterized ellipsoid.
# 
# Roger Stafford 