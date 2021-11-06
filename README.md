# Hairy Black Hole

Este es un paquete que usa numpy, math y scipy para calcular algunos parametros importantes de un agujero negro con pelo estatico. Los parametros que calculo son por ejemplo, la masa, la posicion del horizonte, tanto en x como en r coordenadas. 

Algunas cosas a tomar en cuenta
- No utiliza el sistema internacional de unidades. 
- La distancia se mide en unidades astronomicas, la masa en masas solares, y el tiempo en años. 

Tambien calcula las geodesicas siguiendo un metodo de RK4, las geodesicas tipo tiempo y tipo null.

Para instanciar un agujero negro importamos black_hole y le damos los 3 parametros necesarios que son nu, eta y alpha.

Para isntanciar una geodesica importamos particula_time_like y le damos los 5 parametros necesarios que son nu, eta, alpha, la energia, el momentum angular (todo en las unidades antes mencionadas). En el caso null importamos particula_null y le damos los 4 parametros necesarios que son nu, eta, alpha, b (parametro de impacto).

Esta es una primera version que se ira mejorando, seguramente presentara varios errores con ciertos parametros.