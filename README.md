# Hairy Black Hole

Este es un paquete que usa numpy, math y scipy para calcular algunos parametros importantes de un agujero negro con pelo estático. Los parámetros que calculo por ejemplo son: la masa, la posicion del horizonte, tanto en x como en r coordenadas. 

Algunas cosas a tomar en cuenta
- No utiliza el sistema internacional de unidades. 
- La distancia se mide en unidades astronómicas, la masa en masas solares, y el tiempo en años. 

Tambien calcula las geodésicas siguiendo un metodo de RK4, las geodésicas tipo tiempo y tipo null.

Para instanciar un agujero negro importamos black_hole y le damos los 3 parámetros necesarios que son nu, eta y alpha.

Para instanciar una geodésica importamos particula_time_like y le damos los 5 parametros necesarios que son nu, eta, alpha, la energia, el momentum angular (todo en las unidades antes mencionadas). En el caso null importamos particula_null y le damos los 4 parametros necesarios que son nu, eta, alpha, b (parámetro de impacto).

Esta es una primera versión que se ira mejorando, seguramente presentara varios errores con ciertos parámetros.