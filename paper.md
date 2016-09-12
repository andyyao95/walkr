---
title: 'walkr: MCMC Sampling from Non-Negative Convex Polytopes'
tags:
  - Monte Carlo Markov Chain
  - sampling
  - random walks
  - convex polytope
authors:
 - name: Andy Yu Zhu Yao
   orcid: 0000-0003-3898-8782
   affiliation: Williams College
 - name: David Kane
   affiliation: IQSS, Harvard University
date: 10 September 2016
bibliography: paper.bib
---

# Summary

Consider the intersection of two spaces: the complete solution space 
to Ax = b and the N-simplex. The intersection of these two spaces is 
a non-negative convex polytope. The R package walkr samples from this 
intersection using two Monte-Carlo Markov Chain (MCMC) methods: 
hit-and-run [@kannan] and Dikin walk [@vempala]. Walkr also provide tools to examine sample 
quality [@shinystan]. 

MCMC sampling is of great interest in applied statistics, as it is a common approach to sample
data drawn from a theoretical distribution [@gelman]. In application, walkr will be a powerful tool for estimating
expectations for Bayesian statistics. The walkr package will also be found useful by users who are
interested in generating random weight vectors in high dimensions given specific constraints. 

The real world application to MCMC sampling is vast. In the context of finance, we've had users use 
walkr to generate random portfolios satisfying specific constraints. We've also had users use walkr to sample from 
solution spaces obtained empirically from mass spectrometry analysis of proteins, which can provide 
insight into the biological systems of interest. 
Finally, walkr is one of the first open-sourced softwares to implement the Dikin walk. 

# References
