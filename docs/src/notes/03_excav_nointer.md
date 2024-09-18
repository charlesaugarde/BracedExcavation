# Excavation without interface (robust results)

## Excavation stage 1

To exmaine the influence of constitutive models on the model behaviours, the *Hooke3d*, *LSSO*, *DPconst*, and *DPSSO* models are examined for the first stage of excavation

```@julia
julia> PARAMS = PARAMSDT(CONSTM="Hooke3d")  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav1_hooke3d.png" width="400"/></p>
<p style="text-align: center">Figure 1. The wall deflections due to the first excavation with Hooke3d model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="LSSO")  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav1_LSSO.png" width="400"/></p>
<p style="text-align: center">Figure 2. The wall deflections due to the first excavation with LSSO model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="DPconst", TOL=1e-4, THETA=π/6)  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav1_DPconst_plastic.png" width="400"/></p>
<p style="text-align: center">Figure 3. The wall deflections due to the first excavation with DPconst model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="DPSSO", TOL=1e-4, THETA=π/6)  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav1_DPSSO_plastic.png" width="400"/></p>
<p style="text-align: center">Figure 4. The wall deflections due to the first excavation with DPSSO model.</p>
```
Notes: if the yield function of Drucker-Prager is changed, then the behaviour will change dramatically. 

## Excavation stage 2


```@julia
julia> PARAMS = PARAMSDT(CONSTM="Hooke3d")  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav2_hooke3d.png" width="400"/></p>
<p style="text-align: center">Figure 5. The wall deflections due to the second excavation with Hooke3d model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="LSSO")  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav2_LSSO.png" width="400"/></p>
<p style="text-align: center">Figure 6. The wall deflections due to the second excavation with LSSO model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="DPconst", TOL=1e-4, THETA=π/6)  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav2_DPconst_plastic.png" width="400"/></p>
<p style="text-align: center">Figure 7. The wall deflections due to the second excavation with DPconst model.</p>
```

```@julia
julia> PARAMS = PARAMSDT(CONSTM="DPSSO", TOL=1e-4, THETA=π/6)  
```

```@raw html
<p style="text-align: center"><img src="../figures/nointer_excav2_DPSSO_plastic.png" width="400"/></p>
<p style="text-align: center">Figure 8. The wall deflections due to the second excavation with DPSSO model.</p>
```
