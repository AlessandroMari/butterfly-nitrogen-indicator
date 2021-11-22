# Butterfly nitrogen indicator

### Context and aim
In the context of N deposition, it is crucial to develop an indicator to track the impact of N excess on biodiversity. In this regard, butterflies have many indicator-features as they are well-documented, have a short life cycle, and, hence, quickly react to environmental changes. Moreover, they are deemed to be representative of many terrestrial insect groups. Butterfly optimal nitrogen values, as developed for some Dutch butterfly species in Oostermeijer and Van Swaay (1998), represent a potential candidate as they are based on optimal Ellenberg N indicators for plants, which are closely correlated to the actual N availability in the soil (WallisDeVries & Van Swaay, 2017).

Here, in a three-step R-based process, we first develop a GAM model for calculating the optimal nitrogen value for butterfly species based on those available for Dutch species from Oostermeijer and Van Swaay (1998). The information on butterfly species traits, adopted as principal components from WallisDeVries (2014), is also included in the model implementation. In the second step, the model's predictive ability is then cross-validated on the same Dutch species for which the optimal nitrogen value is available and then comparable. Finally, in the third step, the model is applied to butterfly species occurring in Western European transects to calculate their optimal nitrogen values.

More practical details on how each step is implemented in R will be soon found in the README.md file located in the folder relative to each step.
The complete, but old version of the work as my MSc thesis can be read for free at: https://drive.google.com/file/d/1H_WMqdIaqHgq6hihncuJ07aYCRgjwacu/view?usp=sharing.

#### Additional notes:
- The old version of the project, which used the the principal components from Essens et al. (2017) is also left commented in the script. The new components provide a much neater separation of butterfly traits and therefore have a more significant effect as predictors of the GAM model.
- Example data or documents is not available in the "Data" and "Documents" folder for privacy reasons. Therefore, codes are not reproducible by a github user that accesses this project. However, if data is needed, send an email to "ale_seas@libero.it" and a request to the owner of the data can be delivered. 

#### Citations
- Essens, T., Langevelde, F., Vos, R., Swaay, C., & Wallisdevries, M. (2017). Ecological determinants of butterfly vulnerability across the European continent. Journal of Insect
42
Conservation, 1–12. https://doi.org/10.1007/s10841-017-9972-4
- Oostermeijer, J. G. B., & Van Swaay, C. A. M. (1998). The relationship between butterflies and environmental indicator values: a tool for conservation in a changing landscape. Biological Conservation, 86(3), 271–280.
- WallisDeVries, M. F. (2014). Linking species assemblages to environmental change: Moving beyond the specialist-generalist dichotomy. Basic and Applied Ecology, 15(4), 279–287. https://doi.org/https://doi.org/10.1016/j.baae.2014.05.001
- WallisDeVries, M. F., & Van Swaay, C. A. M. (2017). A nitrogen index to track changes in butterfly species assemblages under nitrogen deposition. Biological Conservation, 212, 448–453.
