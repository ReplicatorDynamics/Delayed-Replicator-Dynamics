clear

Fmax=6
Fstep=100
dstep=80
clf;

F=0:Fmax/Fstep:Fmax
d=0:1/dstep:1
Psi=0.0115
Phi=0.1
alpha=0.0001022    // simplified phenomenological parameter of the juvenile survival term with delay, alpha=-ln(u^(1/z)) 


for z=1:Fstep+1
    for y=1:dstep+1




Delta=(Phi/F(z))^2+Phi/F(z)+0.25-2*Psi/d(y)                      ////(d^2)*(0.5*F+Psi)^2-2*(F^2)*Psi*d
if Delta>=0 then

qstab=0.5-(Phi/F(z))+sqrt(Delta)

Dsup=(((1-qstab)^2)*0.5*d(y)+Psi)/(0.5*F(z)+Phi)       ///((0.5*F+Psi)*d+-sqrt(Delta))/(F^2)           ///!!!!!!!!!!!!!!!!!!!

n=-log(Dsup)/alpha

if n>0 then

n1=qstab*n                      //////(1-(F/d)*Dsup)*n

b11=-0.5*d(y)*(n^2-n1^2)/n^2-Psi;
b12=(0.5*F(z)+Phi)*(1-alpha*n)*exp(-alpha*n);
b13=-d(y)*(1-n1/n);
b14=0;


b21=0;
b22=-(0.5*F(z)*(n1^2)/(n^2)+0.5*F(z)*alpha*(n1^2)/n+alpha*Phi*n1)*exp(-alpha*n);
b23=-Psi;
b24=(F(z)*n1/n+Phi)*exp(-alpha*n);

A=b12*b24-b22*b14
B=b11*b24+b12*b23-b13*b22-b14*b21
C=b11*b23-b13*b21

P=2*((b11+b23)^2)-(b12+b24)^2-4*C
Q=[2*C-(b11+b23)^2]^2+2*(b12+b24)*[(A+C)*(b12+b24)-B*(b11+b23)]+2*(C^2)-2*(A^2)-((b11+b23)*(b12+b24)-B)^2
R=2*(A^2-C^2)*[2*C-(b11+b23)^2]+2*B*(A-C)*[(b11+b23)*(b12+b24)-B]-[(A+C)*(b12+b24)-B*(b11+b23)]^2
S=(A^2-C^2)^2-(B^2)*([A-C]^2)


x = poly(0, "x");
p = S+R*x+Q*x^2+P*x^3+x^4

rootvec=roots(p)

for i=1:4
   if imag(rootvec(i))~=0 then rootvec2(i)=real(0); else rootvec2(i)=real(rootvec(i)) end
end

//u=1

///bifdelay=0

for i=1:4
    u=rootvec2(i)
    if rootvec2(i)>0 
        then
     aux(i)=(1/(sqrt(u)))*acos((-B*(A-C+u)+(b11+b23)*(b12+b24)*u)/(A^2-(C-u)^2-[(b11+b23)^2]*u))
     else aux(i)=0
     end
 end
 
 aux2=real(aux)
 
 maximum=max(aux2)
 
 if maximum<0 then maximum=0 end

for i=1:4
    if aux2(i)<=0 then aux2(i)=maximum end
end

bifdelay(z,y)=min(aux2)

if bifdelay(z,y)<0 then bifdelay(z,y)=0 end

else
    
 bifdelay(z,y)=-0.0001   

end  ////end when n>0

else      //when Delta<0
    
bifdelay(z,y)=--0.0001

end  //end Delta>=0

end    ///end of i iterating
end     ////end of j iterating

height=max(bifdelay)

ebox=[0,1,0,Fmax,0,height]

f1=scf(1);
clf;
plot3d(bifdelay)
//a=gca()
//a.data_bounds=[0,1,0,Fmax,0,height];
//bifdelay=(1/(sqrt(u)))*acos((-B*(A-C)+(b11+b23)*(b12+b24)*u)/(A^2-C^2-[(b11+b23)^2]*u))


