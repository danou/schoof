---
layout: article
title: Preprint
---

Dans ce projet, je vais vous présenter un algorithme de comptage de points de courbes elliptique sur des corps finis. Je me restreindrais à des corps finis $$\mathbf{F}_{p}$$ avec $$p$$ premier différent de 2 et 3. Pour c'est deux derniers cas, l'algorithme est sensiblement le même. 

## Contexte historique

## Courbes elliptiques sur $$\mathbf{F}_{p}$$

Soit $$\mathbf{F}_{p}$$ un corps fini à $$p$$ éléments de caractéristiques $$p\neq 2,3$$.

Soit $$E$$ une courbe elliptique définie sur $$\mathbf{F}_{p}$$. On obtient l'équation affine de Weierstrass : $$y^{2} = x^{3} + ax + b$$ avec $$a,b\in\mathbf{F}_{p}$$ et $$\Delta = -16(4a^{3} + 27b^{2}) \neq 0$$.

1. Definition

Soit $$\varPhi$$ l'endomorphisme de Frobénius d'une courbe elliptique $$E$$ tel que  
$$\begin{array}{clcl}
\varPhi : &E(\bar{\mathbf{F}_{p}}) &\longrightarrow &E(\bar{\mathbf{F}_{p}})\\
&(x, y) &\longmapsto	&(x^{p}, y^{p}).\\
\end{array}$$
