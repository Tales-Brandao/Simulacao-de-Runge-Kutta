import matplotlib.pyplot as plt
import numpy as np
import math


a = [] #concentração celular x
b = [] #tempo 
c = [] #Irradiância 
d = [] #substrato
e = [] #tempo
f = [] #taxa de consumo de enxofre extracelular
g = [] #quota de enxofre no meio extracelular
j = [] #tempo
k = [] #rs
l = [] #taxa de crescimento de celular
m = [] #o2
n = [] #tempo
o = [] #rendimento da produção de o2
p = [] #velocidade de crescimento fotossintetica
r = [] #F(Q)
t = [] #mug 
u = []
v = []

x = 0.07
h = 0.001
D = 0.01
Mus = 0.032
Mup = 0.5
s = 0.0
ys = 11.5
ks = 3.7
mum = 0.2274
s0 = 5.0
yo2 = 1.42
kla = 0.46
o2x = 8.76
o2 = 0.0
ki = 81.38
kii = 2500
q0 = 110
b0 = 0.01728
Ea = 172
Es = 868
fmin = 0.0674
K = 0.3389
qm = 7.0
L = 1
z = 0.5
q = 0

alpha = math.sqrt((Ea/(Ea + 2 * b0 * Es)))
delta = x * (math.sqrt((Ea*(Ea + 2 * b0 * Es))))
rob = math.tanh(delta * (L - z))
rx = (Mup - Mus)* x
rs = ys* mum * ((s*x)/(ks + s))
ro2 = yo2 * rx
rxp = Mup * x
q = (s0 - s)/x
G = (2 * q0 * alpha + rob)/ (((alpha**2)+1)* rob + 2*alpha) 
mug = mum * G / (ki + G + ((G**2)/kii))
fq = fmin + (1 - fmin) * ((math.exp(K-(q/qm)) - 1 )/ (math.exp(K) - 1 ))
print(fq)

def f1(x):
    return ((mug - Mus - D)*x)
 
for i in np.arange (0, 10, 0.001):
    k1 = h*f1(x)
    k2 = h*f1(x+0.5*k1)
    k3 = h*f1(x+0.5*k2)
    k4 = h*f1(x+k3)
    
    x = x + ((k1 + 2*k2 + 2*k3 + k4)/6.0)
    
    delta = x * (math.sqrt((Ea*(Ea + 2 * b0 * Es))))
    rob = math.tanh(delta * (L - z))
    G = (2 * q0 * alpha + rob)/ (((alpha**2)+1)* rob + 2*alpha)
    mug = mum * G / (ki + G + ((G**2)/kii))

    a.append(x)
    b.append(i)
    c.append(G)
    p.append(mug)
    
    #print(delta, G, mug)

def f2(s):
    #print(s)
    return (- ys * mum * ((s*x)/(ks + s)) + D*(s0 - s))

for i in np.arange (0,10, 0.001):
    l1 = h*f2(s)
    l2 = h*f2(s+0.5*l1)
    l3 = h*f2(s+0.5*l2)
    l4 = h*f2(s+l3)
    s = s + ((l1 + 2*l2 + 2*l3 + l4)/6.0)
          
    rs = ys * mum * ((s*x)/(ks + s))
   
    #print(rs,s)

    d.append(s)
    f.append(i)
    e.append(rs)
    
def f3(q):
    return (rs - mug * fq * q)
    
for i in np.arange (0,10, 0.001):
    m1 = h*f3(q)
    m2 = h*f3(q+0.5*m1)
    m3 = h*f3(q+0.5*m2)
    m4 = h*f3(q+m3)
   
    q = q + ((m1 + 2*m2 + 2*m3 + m4)/6.0)
          
    rs = ys * mum * ((s*x)/(ks + s))
    delta = x * (math.sqrt((Ea*(Ea + 2 * b0 * Es))))
    rob = math.tanh(delta * (L - z))
    G = (2 * q0 * alpha + rob)/ (((alpha**2)+1)* rob + 2*alpha)
    mug = mum * G / (ki + G + ((G**2)/kii))
    fq = fmin + (1 - fmin) * ((math.exp(K-(q/qm)) - 1 )/ (math.exp(K) - 1 ))
    
    #print(rs,s, mug)

    g.append(q)
    j.append(i)
    #k.append(rs)
    l.append(rxp)
    r.append(fq)
    t.append(mug)
    #print (t)

def f4(o2):
    return (yo2 * ((mug * fq) - Mus) * x - kla*(o2 - o2x)-D*o2)
    
for i in np.arange (0,10, 0.001):
    n1 = h*f4(o2)
    n2 = h*f4(o2+0.5*n1)
    n3 = h*f4(o2+0.5*n2)
    n4 = h*f4(o2+n3)
   
    o2 = o2 + ((n1 + 2*n2 + 2*n3 + n4)/6.0)
    ro2 = yo2 * rx      
    
    #print(ro2, o2)

    m.append(o2)
    n.append(i)
    o.append(ro2)
    
plt.plot(b,c)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Irradiância")
plt.savefig('Irradiância_por_t.png', format='png')
plt.show() 

plt.plot(b,p)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Velocidade de Crescimento Fotossintética")
plt.savefig('vel_cresc_fotos_por_t.png', format='png')
plt.show() 

plt.plot(b,a)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Concentração de biomassa")
plt.savefig('X_por_t.png', format='png')
plt.show()  

plt.plot(d,a)
#plt.title("Gráfico 1")
plt.xlabel("Concentração de substrato")
plt.ylabel("Concentração de biomassa")
plt.savefig('S_por_X.png', format='png')
plt.show() 

plt.plot(m,a)
#plt.title("Gráfico 1")
plt.xlabel("Concentração de Oxigênio")
plt.ylabel("Concentração de biomassa")
plt.savefig('O2_por_X.png', format='png')
plt.show() 

plt.plot(a,m)
#plt.title("Gráfico 1")
plt.ylabel("Concentração de Oxigênio")
plt.xlabel("Concentração de biomassa")
plt.savefig('X_por_O2.png', format='png')
plt.show()

plt.plot(b,m)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Concentração de Oxigênio")
plt.savefig('O2_por_t.png', format='png')
plt.show() 

plt.plot(b,g)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Quota de enxofre")
plt.savefig('Q_por_t.png', format='png')
plt.show() 

plt.plot(b,l)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Velocidade específica de crescimento celular")
plt.savefig('mi_por_t.png', format='png')
plt.show() 

plt.plot(b,f)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Taxa de consumo de enxofre")
plt.savefig('taxa_enxofre_por_t.png', format='png')
plt.show() 

plt.plot(b,d)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("Substrato")
plt.savefig('S_por_t.png', format='png')
plt.show() 

plt.plot(c,p)
#plt.title("Gráfico 1")
plt.xlabel("Irradiância")
plt.ylabel("Velocidade Específica de Crescimento Celular")
plt.savefig('Irradiância_por_vel_esp_cresc_cel.png', format='png')
plt.show()

plt.plot(g,r)
#plt.title("Gráfico 1")
plt.xlabel("Quota de enxofre")
plt.ylabel("F(Q)")
plt.savefig('Quota_de_enxofre_por_F(Q).png', format='png')
plt.show()

plt.plot(b,r)
#plt.title("Gráfico 1")
plt.xlabel("tempo (dias)")
plt.ylabel("F(Q)")
plt.savefig('Quota_de_enxofre_por_tempo.png', format='png')
plt.show()

plt.plot(c,t)
#plt.title("Gráfico 1")
plt.xlabel("Irradiância")
plt.ylabel("Velocidade de CRescimento na fotossíntese")
plt.savefig('Vel_crssc_fotos_por_Irradiancia.png', format='png')
plt.show()