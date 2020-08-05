# TODO review

- __Linea 574 y 586__: Aplicacion de proporcion de reclutas por sexo, parece errada y sobrestimada respecto de los p_hembras. Verificar curva de proporcion (exp(log_propmR)/(1-exp(log_propmR))).

- __Linea 578 y 591__: Sobreutilizacion de F para el primer ahno. F=1 esta condicionada a la estabilizacion de la estructura para t=1 y la mortalidad por pesca desde los desembarques.

- __Linea 622 y 625__: Buscar consistencia con grupos plus definidos en la condicion inicial, lineas 582 y 594.

- __Eliminar Lineas__: 692, 693, 703. Repite la construccion de la BD


### Task List

:nerd_face:

- [ ] Buscar alcanzar captura 2019 con mayor detalle desde modelo 2018
- [ ] Proyectar con esa captura
- [ ] Graficar con distintas opciones de reclutamiento (opt 1, 2 y 3) observar incidencia de Rec en proyecciones
- [ ] Mismos pasos para LAm_2017
- [ ] Implementar mismo codigo en LAmS.tpl 