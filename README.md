# primitive_calculus

Simply place to check and play with some fun ideas I found

1D nonlinear stuff

Bisection - first order convergence. Terrible method, require specific function and range setup for working

Secant - convergence order is between 1 and 2 (~1.5). Not that bad at all, especially if we don't want take the derivative.

Newton - second order convergence (not always!). Usual dificalies with derivative...

Third order Newton - requires second derivative! It would be much better if we can calculate function derivative analytically. VERY SENSITIVE FOR FUNCTIONS!

It looks like one should solve the problem with optimization method step more carefully in case of ND newton method - simple division by two now frequently does not work.
Yes, if one want to secure convergence, one should solve step minimization problem properly!

The thing with convergence according to normalized defect looks like something rather particular, than general. However this result may be influenced by changes in algorithm of chosing optimal step!

TODO:
[] Try out method of deleting roots for searching all solutions

