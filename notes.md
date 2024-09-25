from https://web.mit.edu/~joshlin/www/dox/deriv__quark__displacement__w_8h_source.html#l00043

https://arxiv.org/pdf/hep-lat/0210030

## derivative quark displacement

$$D_i = s_{ijk}\nabla_j\nabla_k$$
$$B_i = \epsilon_{ijk}\nabla_j\nabla_k$$
the symmetric lattice derivative is:
$$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)$$

where
 
$$s_{ijk} = |\epsilon_{ijk}|$$
$$S_{\alpha jk} = 0\quad j\ne k, S_{111}=S_{222}=1, S_{122}=S_{233}=-1$$
we have set the displacement length(derivative length)to 1 in our chroma input files

in our chroma input files for the meson elementals,which corresponds to the displacement length in derivative
```
disp                     Group
disp_1                   Group
disp_1_1                 Group
disp_1_2                 Group
disp_1_3                 Group
disp_2                   Group
disp_2_1                 Group
disp_2_2                 Group
disp_2_3                 Group
disp_3                   Group
disp_3_1                 Group
disp_3_2                 Group
disp_3_3                 Group
```
* construct the (right Nabla) source 
$\nabla_\mu f(x) = U_\mu(x)f(x+\mu) - U_{-\mu}(x)f(x-\mu)$

$\nabla_\mu F(x)$
    The sink interpolator structure is
    $\Gamma_f \equiv \nabla_i$

* Construct (right D) source

    Operator is  D
    The sink interpolator structure is
    $\Gamma_f \equiv D_i$

* Construct (right B) source

    Operator is  B
    The sink interpolator structure is
    $\Gamma_f \equiv B_i$

* Construct (right E) source

    Operator is  E
    The sink interpolator structure is
    $\Gamma_f \equiv E_i$

* Construct (right Laplacian) source
       Operator is  Laplacian
      The sink interpolator structure is
      $\Gamma_f \equiv Laplacian$
## pion operators 

* Construct (PionxNabla_T1) source
     
    Operator is  Pion x Nabla_T1
    The sink interpolator structure is
    $\Gamma_f \equiv \gamma_5\nabla_i$

* Construct (PionxD_T2) source
     
    Operator is  pion x D_T2
    The sink interpolator is   
    $\Gamma_f \equiv \gamma_4\gamma_5 D_i$  

* Construct (PionxB_T1) source

    Operator is  pion x B_T1
    The sink interpolator is   

    $\Gamma_f \equiv \gamma_5 B_i$  

## mixed gamma operators 

build correlation matrix with multiple gamma operators 

build a 2x2 correlation matrix with both operators gamma_i, nabla_i

derivative based gamma operators result in elementals no longer diagonal in distillation space 

