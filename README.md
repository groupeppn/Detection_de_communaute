# Detection_de_communaute
 L’objectif principal de ce projet est l'implémentation d’algorithmes qui ont pour vocation la détection de communautés dans les réseaux sociaux  représenté en sous forme de graphes , chaque sommets étant un membre et chaque arête une relation entre eux . Dans notre cas  nous avons deux algorithmes à réaliser , ses derniers sont   proposés par notre encadreur  et  font appelle tous deux  à l'algorithme de Bron-Kerbosch . 

 
  Compilation:

 Pour compiler notre programme, on s’est basé sur un Makefile qui regroupe à la fois la partie séquentielle et parallèle.

 Le code en séquentiel est compilé en utilisant la commande make serial et le code parallélisé avec OpenMP est compilé en utilisant la commande make parallel.



Exécution:    

Pour exécuter notre programme  il suffit de taper les commande suivante sur l'invite de commande :

le code séquentiel:    ./codes  -f  mygraph.txt  -a1 


le code parallèle:     ./codep -f  mygraph.txt -t  -a2


Remarque: 

On peut effectuer l'execution en utilisant d'autre graphes qui se trouve dans le dossier graphs.


-f : désigne le nom de fichier avec description du graphe. Les informations sur le graphe sont stockées dans un fichier mygraph1.txt ou nous ajoutons à la fois les arêtes (u, v) et (v, u).

-a indique quel algorithme à exécuter:

 -a 1 : pour l'algorithme 1
 
 -a 2 : pour l’algorithme 2



