function tdot=TemperaturEquationEnsembleTCLs(t,y,s)
Ta=32;
P=1.6e+3;
R=0.0098;
C=74;

tdot(1,1)=(-1/(R*C))*(y - Ta + s*P*R);
end