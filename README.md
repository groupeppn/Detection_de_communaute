# Detection_de_communaute
 L’objectif principal de ce projet est l'implémentation d’algorithmes qui ont pour vocation la détection de communautés dans les réseaux sociaux  représenté en sous forme de graphes , chaque sommets étant un membre et chaque arête une relation entre eux . Dans notre cas  nous avons deux algorithmes à réaliser , ses derniers sont   proposés par notre encadreur  et  font appelle tous deux  à l'algorithme de Bron-Kerbosch . 

 
  Compilation:

 Pour compiler notre programme, on s’est basé sur un Makefile qui regroupe à la fois la partie séquentielle et parallèle.

 Comme le montre la figure ci-dessous, le code en séquentiel est compilé en tapant sur make serial et le code parallélisé avec OpenMP est compilé en tapant sur le terminal make parallel

Exécution    

Pour exécuter notre programme  il suffit de taper les commande suivante sur l'invite de commande :

le code séquentiel:  ./codes -f mygraph.txt -a 1 


le code parallèle:    ./codep -f mygraph.txt - a 2

