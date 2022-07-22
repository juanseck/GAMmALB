%Algoritmo genético para optimizar el Multi-manned Assembly Line Balancing
clear all
clc
%Parámetros de entrada
%Tiempo de ciclo
tc = 138;
%Problema a leer
%nombreArchivo='MITCHELL.IN2';
nombreArchivo='HESKIA.IN2';
%nombreArchivo='MERTENS.IN2';
[n_tareas, t_tareas, n_precedentes, l_precedentes] = leerDatosProblema(nombreArchivo);

%AQUI VA LA PARTE DEL GENETICO
%INICIALIZAR PARAMETROS
%Numero de soluciones
num_sol=25;
%Numero de generaciones
num_gen=100;
%Numero de soluciones elitistas
num_el=3;
%Numero de competidores para torneo
num_comp=2;
%Prob de cruce
prob_cruce=0.75;
%Prob de mutacion
prob_mut=0.25;
%Pesos para el costo del NS, NT, TM
peso_NS=1;
peso_NT=1;
peso_TM=1;

%SOLUCIONES ALEATORIAS
poblacion=zeros(num_sol,n_tareas);
costo=zeros(num_sol,1);

for i=1:num_sol
    poblacion(i,:)=randperm(n_tareas);
end

%EVALUARLAS EN LA FUNCION COSTO
for j=1:num_sol
    %[costo,NS,n_trabajadores,l_estaciones_trabajadores,l_trabajadores_estacion,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,t_muerto] = costoBalanceo(solucion, tc, n_tareas, t_tareas, n_precedentes, l_precedentes,1,1,1);
    costo(j) = costoBalanceo(poblacion(j,:), tc, n_tareas, t_tareas, n_precedentes, l_precedentes,peso_NS,peso_NT,peso_TM);
end

%CICLO DE GENERACIONES
for i=1:num_gen

    %SELECCION
    %Seleccion por elitismo y torneo
    [poblacion, costo] = seleccion(poblacion, costo, num_sol, num_el, num_comp);
    
    %CRUCE
    %Cruce en dos puntos
    [poblacion] = cruce(poblacion, num_sol, n_tareas, prob_cruce, num_el);
    
    %MUTACION
    %Mutación revolviendo la solución con una permutación
    [poblacion] = mutacion(poblacion, num_sol, n_tareas, prob_mut, num_el);
    
    %EVALUARLAS EN LA FUNCION COSTO
    for j=1:num_sol
        %[costo,NS,n_trabajadores,l_estaciones_trabajadores,l_trabajadores_estacion,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,t_muerto] = costoBalanceo(solucion, tc, n_tareas, t_tareas, n_precedentes, l_precedentes,1,1,1);
        costo(j) = costoBalanceo(poblacion(j,:), tc, n_tareas, t_tareas, n_precedentes, l_precedentes,peso_NS,peso_NT,peso_TM);
    end

end

%DESPUES DEL CICLO, ELEGIR AL MEJOR
[mejor_costo,ind_mejor]=min(costo);

%Imprimir solucion
[costo,NS,n_trabajadores,l_estaciones_trabajadores,l_trabajadores_estacion,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,t_muerto,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt] = costoBalanceo(poblacion(ind_mejor,:),tc,n_tareas,t_tareas,n_precedentes,l_precedentes,peso_NS,peso_NT,peso_TM);
disp(['Solucion: ' mat2str(poblacion(ind_mejor,:))]);
disp(['Costo: ' num2str(mejor_costo)]);
disp(['Numero de estaciones: ' num2str(NS)]);
disp(['Numero de trabajadores: ' num2str(n_trabajadores)]);
%Para tener tabla de Gantt final, ordenar primero por trabajador y luego
%por tiempo de inicio
Tabla=sortrows(tabla_Gantt,[3]);
%disp(Tabla)
diagrama_Gantt_trabajadores(Tabla)

%Cromosoma de prueba
%solucion=[1,4,2,5,7,3,6];
%solucion=[6,4,1,2,5,7,3];
%[costo,NS,n_trabajadores,l_estaciones_trabajadores,l_trabajadores_estacion,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,t_muerto] = costoBalanceo(solucion, tc, n_tareas, t_tareas, n_precedentes, l_precedentes,1,1,1);
%costo

function [nuevaPoblacion, nuevoCosto] = seleccion(Poblacion, costo, num_sol, num_el, num_comp)
% Seleccion completa con elitismo y torneo
[nuevaPoblacion, nuevoCosto] = seleccionElitista(Poblacion, costo, num_el);
ind=num_el+1;
while (ind<=num_sol)
    indGanador = seleccionTorneo(costo, num_sol, num_comp);
    nuevaPoblacion(ind,:)=Poblacion(indGanador,:);
    nuevoCosto(ind)=costo(indGanador);
    ind=ind+1;
end

end

function [nuevaPoblacion, nuevoCosto] = seleccionElitista(Poblacion, costo, num_el)
%Selección elitista de numero individuos de Poblacion
nuevaPoblacion=Poblacion;
nuevoCosto=costo;
[ ~, indices ] = sort(costo);
for i=1:num_el
    nuevaPoblacion(i,:)=Poblacion(indices(i),:);
    nuevoCosto(i)=costo(indices(i));
end
end

function indGanador= seleccionTorneo(costo,num_sol,num_comp)
%Selecciona un ganador por torneo entre numCompetidores aleatorios de Poblacion
indices=randperm(num_sol,num_comp);
[~,ind]=min(costo(indices));
indGanador=indices(ind);
end

function [Poblacion] = cruce(Poblacion, num_sol, n_tareas, prob_cruce, num_el)

% Cruce HA, para operaciones 50% POX y 50% JBX, para maquinas cruce Dos
% Puntos, respeta las primeras numElitista soluciones por elitismo
% Lista de cruces aleatorios
if mod(num_sol-num_el,2)==1
    num_ef=num_el+1;
else
    num_ef=num_el;
end
listaC=randperm(num_sol-num_ef)+num_ef;
listaC=reshape(listaC,(num_sol-num_ef)/2,2);
for i=1:(num_sol-num_ef)/2
    probC=rand;
    if probC<=prob_cruce
        solucion1_so=Poblacion(listaC(i,1),:);
        solucion2_so=Poblacion(listaC(i,2),:);
        %Cruce operaciones
        %probCO=rand;
        %if probCO<=0.5
        %    [solucion1_so,solucion2_so] = crucePOX(num_trab, solucion1_so, solucion2_so);
        %else
        %    [solucion1_so,solucion2_so] = cruceJBX(num_trab, solucion1_so, solucion2_so);
        %end
        [solucion1_so,solucion2_so] = cruceJBX(n_tareas,solucion1_so,solucion2_so);
        Poblacion(listaC(i,1),:)=solucion1_so;
        Poblacion(listaC(i,2),:)=solucion2_so;
    end
end
end

function [nuevaSol1,nuevaSol2] = crucePOX(numtrabajos, solucion1, solucion2)
%Cruce POX de secuencias de operaciones
%Se divide el conjunto de trabajos en dos grupos g1 y g2 de forma aleatoria
%Cualquier elemento de secuencia1 que pertenezca a g1 se acomoda en orden1 en las mismas posiciones
%el resto de las posiciones se llenan con los elementos de g2 en el orden en que aparecen en secuencia2
%De manera similar se construye orden2 interca,biando los papeles de secuencia1 y secuencia2
nuevaSol1= solucion1;
nuevaSol2= solucion2;
%Numero de trabajos en g1, se genera numero aleatorio entre 1 y
%numtrabajos-1
numg1=randi([1,numtrabajos-1]);
%Permutacion aleatoria de trabajos
per=randperm(numtrabajos);
g1=per(1:numg1);
g2=per(numg1+1:numtrabajos);
[~,indices1]=ismember(solucion1,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion2,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol1(indices1N)=solucion2(indices2);
%Para la solucion 2
[~,indices1]=ismember(solucion2,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion1,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol2(indices1N)=solucion1(indices2);
end

function [nuevaSol1,nuevaSol2] = cruceJBX(numtrabajos, solucion1, solucion2)
%Cruce JBX de secuencias de operaciones
%Se divide el conjunto de trabajos en dos grupos g1 y g2 de forma aleatoria
%Cualquier elemento de secuencia1 que pertenezca a g1 se acomoda en orden1 en las mismas posiciones
%el resto de las posiciones se llenan con los elementos de g2 en el orden en que aparecen en secuencia2
%Cualquier elemento de secuencia2 que pertenezca a g2 se acomoda en orden2 en las mismas posiciones,
%el resto de las posiciones se llenan con los elementos de g1 en el orden en que aparecen en secuencia1
nuevaSol1= solucion1;
nuevaSol2= solucion2;
%Numero de trabajos en g1, se genera numero aleatorio entre 1 y
%numtrabajos-1
numg1=randi([1,numtrabajos-1]);
%Permutacion aleatoria de trabajos
per=randperm(numtrabajos);
g1=per(1:numg1);
g2=per(numg1+1:numtrabajos);
[~,indices1]=ismember(solucion1,g1);
indices1=indices1>0;
[~,indices2]=ismember(solucion2,g2);
indices2=indices2>0;
indices1N=(~indices1);
nuevaSol1(indices1N)=solucion2(indices2);
%Para la solucion 2
indices2N=(~indices2);
nuevaSol2(indices2N)=solucion1(indices1);
end

function [nuevaSol1,nuevaSol2] = cruceDosPuntos(numOperaciones,solucion1,solucion2)
%Esta función calcula el cruce en dos puntos aleatorios de dos secuencias de máquinas
%Se divide cada secuencia en dos puntos aleatorios intermedios
%y se intercambian las secuencias centrales para formar secuencias nuevas
%Para el algoritmo HA, como las secuencias de máquinas ya son factibles, y no se
%cambian posiciones entre las soluciones padres, entonces las secuencias generadas
%también serán factibles, por lo que no hay necesidad de una revisión extra.
nuevaSol1= solucion1;
nuevaSol2= solucion2;
%Lugares aleatorios de corte
per=randperm(numOperaciones);
lugar1=per(1);
if lugar1>per(2)
    lugar1=per(2);
    lugar2=per(1);
else
    lugar2=per(2);
end
%Hacer cruce
nuevaSol1(1:lugar1)=solucion1(1:lugar1);
nuevaSol1(lugar1+1:lugar2)=solucion2(lugar1+1:lugar2);
nuevaSol1(lugar2+1:numOperaciones)=solucion1(lugar2+1:numOperaciones);
nuevaSol2(1:lugar1)=solucion2(1:lugar1);
nuevaSol2(lugar1+1:lugar2)=solucion1(lugar1+1:lugar2);
nuevaSol2(lugar2+1:numOperaciones)=solucion2(lugar2+1:numOperaciones);
end

function [Poblacion] = mutacion(Poblacion,num_sol,n_tareas,prob_mut,num_el)
% Cruce HA, para operaciones 50% Intercambio y 50% Vecindad, para maquinas
% mutacionMaquinas
for i=num_el+1:num_sol
    probM=rand;
    if probM<=prob_mut
        %Mutación operaciones
        Poblacion(i,:) = mutacionIntercambio(n_tareas,Poblacion(i,:));
        %probMO=rand;
        %if probMO<=0.5
        %    Poblacion(i,:) = mutacionIntercambio(n_tareas,Poblacion(i,:));
        %else
        %    Poblacion(i,:) = mutacionVecindad(numeroTrabajos,n_tareas,Poblacion(i,:));
        %end
    end
end

end

function [nuevaSol] = mutacionIntercambio(numOperaciones,solucion)
%Mutación de operaciones
%Es muy sencilla, lo que hace es intercambiar solo dos posiciones aleatorias de secuencia de operaciones
%para generar la nueva orden de operaciones
nuevaSol = solucion;
%Posiciones a intercambiar
pos=randperm(numOperaciones,2);
nuevaSol(pos(1))=solucion(pos(2));
nuevaSol(pos(2))=solucion(pos(1));
end

function [nuevaSol] = mutacionVecindad(numTrabajos,numeroOperaciones,solucion)
%Mutación de operaciones
%Consiste en seleccionar 3 elementos aleatorios de secuencia, los elementos deben ser distintos, de distintos trabajos.
%para generar la nueva orden de operaciones
%Hacer todos los posibles intercambios de posiciones de esos 3 elementos (vecinos) y seleccionar uno de esos vecinos al azar como nuevo orden
nuevaSol = solucion;
posiciones=1:numeroOperaciones;
posicionesAl=zeros(1,3);
%Tres aleatorios trabajos
trab=randperm(numTrabajos,3);
%Tomar posiciones aleatorias de cada trabajo
for i=1:3
    indices=solucion==trab(i);
    pos=posiciones(indices);
    posicionesAl(i)=pos(randi(length(pos),1));
end
%Permutar posiciones
aux=randperm(3);
posicionesFin=posicionesAl(aux);
%Hacer vecino
nuevaSol(posicionesFin)=solucion(posicionesAl);
end

function [costo,NS,n_trabajadores,l_estaciones_trabajadores,l_trabajadores_estacion,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,t_muerto,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt] = costoBalanceo(solucion, tc, n_tareas, t_tareas, n_precedentes, l_precedentes,factorNS,factorNT,factorTM)
%Inicializar costo del balanceo
%costo=0;
%Mínimo número de trabajadores por estación
w_min=1;
%Máximo número de trabajadores por estación
w_max=2;
%Parametro de control de cota superior
beta=2;
%Parámetro de castigo por no cumplir la precedencia de tareas
M=1e6;

%Algoritmo
%Calcular Wlb
Wlb=round(sum(t_tareas)/tc);
%Calcular UB
UB=beta*((tc*Wlb-sum(t_tareas))/Wlb);
%Número inicial de estaciones de trabajo
NS=0;
%Arreglos para el control de trabajadores en estaciones y tiempos de cada trabajador
n_trabajadores=0;
%Lista de trabajadores en cada estacion
l_estaciones_trabajadores=zeros(n_tareas,n_tareas);
%Lista de numero de trabajadores por estación
l_trabajadores_estacion=zeros(n_tareas,1);
%Tiempo de trabajo acumulado de cada trabajador
t_trabajador=zeros(n_tareas,1);
%Lista del tiempo de finalización de cada trabajador
l_tiempofin_trabajador=zeros(n_tareas,1);
%Lista de tiempo de finalización de cada tarea
l_tiempofin_tarea=zeros(n_tareas,1);
%Lista de trabajos asignados a cada trabajador
l_trabajadores_trabajos=zeros(n_tareas,n_tareas);
%Numero de trabajos asignados a cada trabajador
n_trabajadores_trabajos=zeros(n_tareas,1);
%Tabla Gantt, los renglones son tareas, las columnas son (tarea,
%trabajador, inicio, duracion)
tabla_Gantt=zeros(n_tareas,5);
%Cromosoma principal
c_pri=solucion;
%Cromosoma auxiliar que guarda tareas que no han cumplido su precedencia
n_aux=0;
c_aux=zeros(1,n_tareas);
%Tiempo muerto
t_muerto=zeros(n_tareas,1);

%HACER DOS CICLOS FOR, UNO PARA EL CROMOSOMA PRINCIPAL Y OTRO PÀRA LAS
%TAREAS QUE NO SE HAYAN PODIDO ACOMODAR POR PRECEDENCIA


%Recorrer el cromosoma
for i=1:n_tareas
    %Tarea por revisar en el cromosoma principal
    tarea=c_pri(i);
    %Tiempo de la tarea
    tt=t_tareas(tarea);
    %Funcion para revisar si hay tareas precedentes no realizadas o no
    [bandera_precedencia,FTP]=tareas_precedentes(tarea,n_precedentes,l_precedentes,l_tiempofin_tarea);
    %Si la tarea cumple la precedencia, tratar de acomodarla
    if bandera_precedencia == 1
        %Acomodar la tarea en estaciones actuales o una nueva
        [NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt]=acomodar_tarea(tarea,w_max,tc,tt,FTP,NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt);
    else %No todas las tareas precedentes se han hecho, poner la tarea en un cromosoma auxiliar
        n_aux=n_aux+1;
        c_aux(n_aux)=tarea;
    end
end

%Si hay tareas que no se pudieron acomodar por precedencia, acomodarlas
%ahora
while n_aux >0
    %Encontrar tarea que ya se pueda acomodar
    for i=1:n_aux
        %Tarea a revisar
        tarea=c_aux(i);
        %Tiempo de la tarea
        tt=t_tareas(tarea);
        %Funcion para revisar si hay tareas precedentes no realizadas o no
        [bandera_precedencia,FTP]=tareas_precedentes(tarea,n_precedentes,l_precedentes,l_tiempofin_tarea);
        %Si la tarea cumple la precedencia, tratar de acomodarla
        if bandera_precedencia == 1
            %Acomodar la tarea en estaciones actuales o una nueva
            [NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt]=acomodar_tarea(tarea,w_max,tc,tt,FTP,NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt);
            %Quitar tarea de c_aux, decrementar n_aux y salir del ciclo
            c_aux(i)=[];
            n_aux=n_aux-1;
            break;
        end
    end
end

%Calcular tiempo muerto de cada estación
for i=1:NS
    %Trabajadores por estacion
    trabajadores=l_estaciones_trabajadores(i,1:l_trabajadores_estacion(i));
    tiempo_ocupado=sum(t_trabajador(trabajadores));
    t_muerto(i)=((l_trabajadores_estacion(i)*tc)-tiempo_ocupado)/tc;
end
%Calular el costo del balanceo 
costo = (factorNS*NS)+(factorNT*n_trabajadores);
for i=1:NS
    if l_trabajadores_estacion(i)>w_min && t_muerto(i)>UB
        costo = costo + (t_muerto(i)-UB)*factorTM;
    end
end
end
 
function [n_tareas, t_tareas, n_precedentes, l_precedentes] = leerDatosProblema(nombreArchivo)
%La función lee la tabla de datos del problema 
%   Leer archivo
archivo=fopen(nombreArchivo,'r');
%   Buscar todos los datos numéricos enteros que definen el problema
datos=fscanf(archivo,'%c');
% Leer datos hasta encontrar un espacio
in_i=1;
in_f=in_i;
while(datos(in_f) ~= 13)
        in_f=in_f+1;
end
n_tareas=str2num(datos(in_i:in_f));
% Obtener datos tiempos de tareas
t_tareas=zeros(1,n_tareas);
for i=1:n_tareas
    in_i=in_f+1;
    aux=datos(in_i);
    while(aux == 10)
        in_i=in_i+1;
        aux=datos(in_i);
    end
    in_f=in_i;
    aux=datos(in_f);
    while(aux ~= 13)
        in_f=in_f+1;
        aux=datos(in_f);
    end
    t_tareas(i)=str2num(datos(in_i:in_f-1));
end
% Inicializar vectores de número y lista de tareas precedentes
n_precedentes=zeros(1,n_tareas);
l_precedentes=zeros(n_tareas,n_tareas);
% Empezar a leer tareas precedentes
while(1)
    in_i=in_f+1;
    aux=datos(in_i);
    %Tarea precedente
    while(aux == 10)
        in_i=in_i+1;
        aux=datos(in_i);
        if aux=="-"
            break
        end
    end
    if aux=="-"
        break
    end
    in_f=in_i;
    aux=datos(in_f);
    while(aux ~= 44)
        in_f=in_f+1;
        aux=datos(in_f);
    end
    tarea_i=str2num(datos(in_i:in_f-1));
    %Tarea siguiente
    in_i=in_f+1;
    aux=datos(in_i);
    while(aux == 10)
        in_i=in_i+1;
        aux=datos(in_i);
    end
    in_f=in_i;
    aux=datos(in_f);
    while(aux ~= 13)
        in_f=in_f+1;
        aux=datos(in_f);
    end
    tarea_f=str2num(datos(in_i:in_f-1));
    %Tomar las tareas precedentes de tarea_i
    n_aux=n_precedentes(tarea_i);
    l_aux=l_precedentes(tarea_i,1:n_aux);
    %Agregarlas como tareas precedentes de tarea_f si no están
    for i=1:n_aux
        bandera_tarea=0;
        for j=1:n_precedentes(tarea_f)
            if l_precedentes(tarea_f,j)==l_aux(i)
                bandera_tarea=1;
                break
            end
        end
        if bandera_tarea==0
            n_precedentes(tarea_f)=n_precedentes(tarea_f)+1;
            l_precedentes(tarea_f,n_precedentes(tarea_f))=l_aux(i);
        end
    end
    %Agregar tarea_i como precedente de tarea_f
    n_precedentes(tarea_f)=n_precedentes(tarea_f)+1;
    l_precedentes(tarea_f,n_precedentes(tarea_f))=tarea_i;
end
%Maximo de número de tareas precedentes
max_nt=max(n_precedentes);
%Cortar las columnas de l_precedentes no usadas
l_precedentes(:,max_nt+1:end)=[];
end

%Funcion para revisar si hay tareas precedentes no realizadas o no
function [bandera_precedencia,FTP]=tareas_precedentes(tarea,n_precedentes,l_precedentes,l_tiempofin_tarea)
%Bandera que indica si la tarea seleccioanda cumple la precedencia
bandera_precedencia=0;
FTP=0;
%Numero de tareas precedentes de tarea
aux1=n_precedentes(tarea);
if aux1==0 %No hay tareas precedentes, el tiempo de la tarea precedente es cero
    bandera_precedencia=1;
else
    %Revisar si ya se realizaron las actividades precedentes
    %Lista de tareas precedentes
    aux2=l_precedentes(tarea,1:aux1);
    %Ver si alguna tarea precedente tiene tiempo cero, entonces no se
    %ha programado y la tarea actual no se puede hacer
    aux3=l_tiempofin_tarea(aux2);
    if sum(aux3==0)==0 %Las tareas precedentes ya se hicieron, poner el tiempo mayor de estas tareas
        FTP=max(aux3);
        bandera_precedencia=1;
    end
end
end

%Funcion para acomodar la tarea en estaciones de trabajo
function [NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt]=acomodar_tarea(tarea,w_max,tc,tt,FTP,NS,n_trabajadores,l_trabajadores_estacion,l_estaciones_trabajadores,t_trabajador,l_tiempofin_trabajador,l_tiempofin_tarea,l_trabajadores_trabajos,n_trabajadores_trabajos,tabla_Gantt)
%Bandera que indica si la tarea seleccioanda se pudo acomodar en una estación
bandera_tarea=0;
%Recorrer estaciones ya creadas para ver en cual se puede acomodar la tarea
for j=1:NS
    %Revisar los trabajadores de la estación j
    for k=1:l_trabajadores_estacion(j)
        %Trabajador k en la estación j
        trabajador=l_estaciones_trabajadores(j,k);
        %Tiempo final del trabajador k
        TWL=t_trabajador(trabajador);
        %Revisar si se respeta el tiempo de ciclo
        %if TWL+tt <= tc && FTP+tt <= tc %No se usa el FTP pero se deja como referencia
        if TWL+tt <= tc
            tiempo_final_trabajador=TWL+tt;
            tiempo_inicial_tarea=max(l_tiempofin_trabajador(trabajador),FTP);
            tiempo_final_tarea=tiempo_inicial_tarea+tt;
            trabajador_final=trabajador;
            estacion_final=j;
            %Revisar si no había otro trabajador ya seleccionado
            if bandera_tarea==0
                bandera_tarea=1;
                tiempo_seleccionado_trabajador=tiempo_final_trabajador;
                tiempo_seleccionado_inicio_tarea=tiempo_inicial_tarea;
                tiempo_seleccionado_tarea=tiempo_final_tarea;
                trabajador_seleccionado=trabajador_final;
                estacion_seleccionada=estacion_final;
            else
                %Ver si el tiempo final es menor al ya seleccionado
                if tiempo_final_tarea<tiempo_seleccionado_tarea
                    tiempo_seleccionado_trabajador=tiempo_final_trabajador;
                    tiempo_seleccionado_inicio_tarea=tiempo_inicial_tarea;
                    tiempo_seleccionado_tarea=tiempo_final_tarea;
                    trabajador_seleccionado=trabajador_final;
                    estacion_seleccionada=estacion_final;
                elseif  tiempo_final_tarea == tiempo_seleccionado_tarea
                    %Tirar un volado para ver si se sustituye el
                    %trabajador y la estación seleccionadas
                    if rand <0.5
                        tiempo_seleccionado_trabajador=tiempo_final_trabajador;
                        tiempo_seleccionado_inicio_tarea=tiempo_inicial_tarea;
                        trabajador_seleccionado=trabajador_final;
                        estacion_seleccionada=estacion_final;
                    end
                end
            end
        end
    end
end
if bandera_tarea==1 %Si se acomoda la tarea, actualizar tiempos finales de trabajador
    t_trabajador(trabajador_seleccionado)=tiempo_seleccionado_trabajador;
    l_tiempofin_trabajador(trabajador_seleccionado)=tiempo_final_tarea;
    l_tiempofin_tarea(tarea)=tiempo_final_tarea;
    n_trabajadores_trabajos(trabajador_seleccionado)=n_trabajadores_trabajos(trabajador_seleccionado)+1;
    l_trabajadores_trabajos(trabajador_seleccionado,n_trabajadores_trabajos(trabajador_seleccionado))=tarea;
    tabla_Gantt(tarea,1)=tarea;
    tabla_Gantt(tarea,2)=trabajador_seleccionado;
    tabla_Gantt(tarea,3)=tiempo_seleccionado_inicio_tarea;
    tabla_Gantt(tarea,4)=tt;
    tabla_Gantt(tarea,5)=estacion_seleccionada;
else %Si no se puedo acomodar la tarea, se asigna un nuevo trabajador a la última estacióm o se abre una estación nueva
    if NS==0 || l_trabajadores_estacion(NS)==w_max
        NS=NS+1;
    end
    %poner nuevo trabajador en la estación
    n_trabajadores=n_trabajadores+1;
    l_trabajadores_estacion(NS)=l_trabajadores_estacion(NS)+1;
    l_estaciones_trabajadores(NS,l_trabajadores_estacion(NS))=n_trabajadores;
    t_trabajador(n_trabajadores)=tt;
    l_tiempofin_trabajador(n_trabajadores)=FTP+tt;
    l_tiempofin_tarea(tarea)=FTP+tt;
    n_trabajadores_trabajos(n_trabajadores)=n_trabajadores_trabajos(n_trabajadores)+1;
    l_trabajadores_trabajos(n_trabajadores,n_trabajadores_trabajos(n_trabajadores))=tarea;
    tabla_Gantt(tarea,1)=tarea;
    tabla_Gantt(tarea,2)=n_trabajadores;
    tabla_Gantt(tarea,3)=FTP;
    tabla_Gantt(tarea,4)=tt;
    tabla_Gantt(tarea,5)=NS;
end
end

%Diagrama de gantt de los trabajadores
function diagrama_Gantt_trabajadores(Tabla)

%Grafica el diagrama de Gantt de una solucion para el problema de FJSP
%maxNumOper=max(operActualMaquina);

%Numero de trabajos
nt=size(Tabla,1);

%Colores del diagrama
%colores_trab=rand(nt,3);
num_trab=max(Tabla(:,2));
colores_trab = linspecer(nt);
%colores_trab = linspecer(nt);
%colores=(0:(1/num_trab):1)+rand;
%colores(colores>1)=colores(colores>1)-1;
% colores=rand(1,num_trab);
% colores_trab=colores';
% for i=1:2
%     colores2=colores+rand;
%     colores2(colores2>1)=colores2(colores2>1)-1;
%     colores_trab=[colores_trab colores2'];
% end
%colorT=colormap(hsv(nt));
%colorT=colormap(lines(nt+1));


num_est=max(Tabla(:,5));
%colores_est=rand(nt,3);
colores_est = linspecer(nt);
% colores=(0:(1/num_est):1)+rand;
% colores(colores>1)=colores(colores>1)-1;
% colores_est=colores';
% for i=1:2
%     colores2=colores+rand;
%     colores2(colores2>1)=colores2(colores2>1)-1;
%     colores_est=[colores_est colores2'];
% end
%colorT=colormap(hsv(nt));
%colorT=colormap(lines(nt+1));


%Grafica del diagrama de trabajadores
figure(1);
colorT=colormap(colores_trab);
clf
grid on
grid minor
%Ciclo para maquinas
for trabajo=1:nt
    %y=[trabajo-0.75 trabajo-0.25 trabajo-0.25 trabajo-0.75];
    %x=[-4 -4 -1 -1 ];
    %patch(x,y,colorM(trabajo,:));
    %etiqueta=strcat( 'M', num2str(trabajo));
    %text(-3.5,trabajo-0.5,etiqueta)
    y=[trabajo-0.9 trabajo-0.1 trabajo-0.1 trabajo-0.9];
    x=[Tabla(trabajo,3) Tabla(trabajo,3) Tabla(trabajo,3)+Tabla(trabajo,4) Tabla(trabajo,3)+Tabla(trabajo,4)];
    patch(x,y,colorT(Tabla(trabajo,2)));
    etiqueta=strcat('W:', num2str(Tabla(trabajo,2)), ', T:', num2str(Tabla(trabajo,1)));
    text(Tabla(trabajo,3),trabajo-0.5,etiqueta)
    %Para cada máquina, recorrer sus operaciones
    %for oper=1:operActualMaquina(trabajo)
    %    x = [tiemposInMaq(trabajo,oper) tiemposInMaq(trabajo,oper) tiemposFinMaq(trabajo,oper) tiemposFinMaq(trabajo,oper)];
    %    patch(x,y,colorT(trabMaq(trabajo,oper),:));
    %    %etiqueta=strcat( num2str(solucion.tablasMaq.trabajo(maquina,oper)), ',', num2str(solucion.tablasMaq.operacion(maquina,oper)), ',', num2str(solucion.tablasMaq.tiempoIn(maquina,oper)), ',', num2str(solucion.tablasMaq.tiempoFin(maquina,oper)) ); 
    %    etiqueta=strcat( 'J', num2str(trabMaq(trabajo,oper)), ',', num2str(operTrabMaq(trabajo,oper)) ); 
    %    text(tiemposInMaq(trabajo,oper),trabajo-0.5,etiqueta)
    %end
end
%Limite eje x
xlim([-4,Tabla(nt,3)+Tabla(nt,4)+5])
%Remover etiquetas en eje y y empezar eje x en 0
grafica=gca;
grafica.YDir = 'reverse';
set(grafica,'YTick',[])
%etiquetasY=get(grafica,'XtickLabel');
%etiquetasY{1}='';
%set(grafica,'XtickLabel',etiquetasY);
%Linea de 0 a numeroMaquinas en eje y para x en 0
x=[0,0];
y=[0,nt];
hold on
plot(x,y,'k','LineWidth',1.25)
title('Diagrama de Gantt de trabjadores')
hold off

%Grafica del diagrama de estaciones
figure(2);
colorE=colormap(colores_est);
clf
grid on
grid minor
%Ciclo para maquinas
for trabajo=1:nt
    %y=[trabajo-0.75 trabajo-0.25 trabajo-0.25 trabajo-0.75];
    %x=[-4 -4 -1 -1 ];
    %patch(x,y,colorM(trabajo,:));
    %etiqueta=strcat( 'M', num2str(trabajo));
    %text(-3.5,trabajo-0.5,etiqueta)
    y=[trabajo-0.9 trabajo-0.1 trabajo-0.1 trabajo-0.9];
    x=[Tabla(trabajo,3) Tabla(trabajo,3) Tabla(trabajo,3)+Tabla(trabajo,4) Tabla(trabajo,3)+Tabla(trabajo,4)];
    patch(x,y,colorE(Tabla(trabajo,5)));
    etiqueta=strcat('S:', num2str(Tabla(trabajo,5)), ', T:', num2str(Tabla(trabajo,1)));
    text(Tabla(trabajo,3),trabajo-0.5,etiqueta)
    %Para cada máquina, recorrer sus operaciones
    %for oper=1:operActualMaquina(trabajo)
    %    x = [tiemposInMaq(trabajo,oper) tiemposInMaq(trabajo,oper) tiemposFinMaq(trabajo,oper) tiemposFinMaq(trabajo,oper)];
    %    patch(x,y,colorT(trabMaq(trabajo,oper),:));
    %    %etiqueta=strcat( num2str(solucion.tablasMaq.trabajo(maquina,oper)), ',', num2str(solucion.tablasMaq.operacion(maquina,oper)), ',', num2str(solucion.tablasMaq.tiempoIn(maquina,oper)), ',', num2str(solucion.tablasMaq.tiempoFin(maquina,oper)) ); 
    %    etiqueta=strcat( 'J', num2str(trabMaq(trabajo,oper)), ',', num2str(operTrabMaq(trabajo,oper)) ); 
    %    text(tiemposInMaq(trabajo,oper),trabajo-0.5,etiqueta)
    %end
end
%Limite eje x
xlim([-4,Tabla(nt,3)+Tabla(nt,4)+5])
%Remover etiquetas en eje y y empezar eje x en 0
grafica=gca;
grafica.YDir = 'reverse';
set(grafica,'YTick',[])
%etiquetasY=get(grafica,'XtickLabel');
%etiquetasY{1}='';
%set(grafica,'XtickLabel',etiquetasY);
%Linea de 0 a numeroMaquinas en eje y para x en 0
x=[0,0];
y=[0,nt];
hold on
plot(x,y,'k','LineWidth',1.25)
title('Diagrama de Gantt de estaciones')
hold off




end