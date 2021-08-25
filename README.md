# Butterfly nitrogen indicator

### Context and aim
In the context of N deposition, it is crucial to develop an indicator to track the impact of N excess on biodiversity. In this regard, butterflies have many indicator-features as they are well-documented, have a short life cycle, and, hence, quickly react to environmental changes. Moreover, they are deemed to be representative of many terrestrial insect groups. Butterfly optimal nitrogen values, as developed for some Dutch butterfly species in Oostermeijer and Van Swaay (1998), represent a potential candidate as they are based on optimal Ellenberg N indicators for plants, which are closely correlated to the actual N availability (WallisDeVries & Van Swaay, 2017).

Here, in a three-step R-based process, we first develop a GAM model for calculating the optimal nitrogen value for butterfly species based on those available for Dutch species from Oostermeijer and Van Swaay (1998). The information on butterfly species traits and habitat specialization, adopted from Essens et al. (2017), is also included in the model implementation. In the second step, the model's predictive ability is then cross-validated on the same Dutch species for which the optimal nitrogen value is available and then comparable. Finally, the model is applied to butterfly species occurring in Western European transects to calculate their optimal nitrogen values.

More practical details on how each step is implemented in R can be found in the README.md file located in the folder relative to each step.

#### Citations
- Essens, T., Langevelde, F., Vos, R., Swaay, C., & Wallisdevries, M. (2017). Ecological determinants of butterfly vulnerability across the European continent. Journal of Insect
42
Conservation, 1–12. https://doi.org/10.1007/s10841-017-9972-4
- Oostermeijer, J. G. B., & Van Swaay, C. A. M. (1998). The relationship between butterflies and environmental indicator values: a tool for conservation in a changing landscape. Biological Conservation, 86(3), 271–280.
- WallisDeVries, M. F., & Van Swaay, C. A. M. (2017). A nitrogen index to track changes in butterfly species assemblages under nitrogen deposition. Biological Conservation, 212, 448–453.