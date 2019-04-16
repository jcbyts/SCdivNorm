# SCdivNorm
Can SC has divisive normalization?


## Basic approach
The goal of this project is to establish whether SC is well described by divisive normalization when multiple targets are presented (without saccade planning) AND to map the divisive field.

We approach the problem probabilistically by defining a generative model of the spiking process in SC and then we try to fit that model.

Spike counts at time bin $t$ are generated from a poisson process with rate $\lambda(t)$ and bin size $\Delta$,

$r(t)\sim Poiss(\lambda(t)\Delta)$, and

$\lambda(t) = f\big[ g_e[k_e\cdot s(t)] \times g_s[k_s \cdot s(t)]\big]$

where $k_e$ is the neuron's receptive field and $k_s$ is the divisive field that we are interested in. $f$ is the spike nonlinearity and $g_e$ and $g_s$ are subunit nonlinearities for the receptive field and suppresive field, respectively.

The receptive field drive, $g_e[k_e \cdot s(t)]$, is multiplied by the divisive drive $g_s[k_s \cdot s(t)]$. Multiplied? Yes, multipled. We will constrain $g_s[.]$ to span $[{0,1}]$ such that it will always be suppressive when multipled by the receptive field drive.

This model is not convex, so there is not one global optimum. This will probably create some problems. We need to make sure the model will converge to some sort of global optimum. The plan is to optimize using coordinate ascent, iteratively fitting $k_e$ and $k_s$. But, we probably have to initialize the filters in a clever way.

Of course, we need some comparison points (models that fail). For these we will also fit
1) Linear Nonlinear Poisson (LNP)
2) Nonlinear Input Model (NIM)
3) Spike Triggered Covariance (STC)

We compare those to divisive suppression (divS).


In this repo, you will find a sub-directory called `analyses/2019_initial_simulations' that contains the file `scnormsim.m`. Edit that file and save a copy so you can play with stuff. 

It lets you tweak the parameters of a model to generate normalization-like behavior from an SC neuron.

The first figure shows the output of the model neuron. 

![alt text][thedata]

You can see that the "neuron" behaves differently depending on the number of targets.

If we fit a linear nonlinear poisson (LNP) model to these data, it does pretty terribly. Although it does get the right location of the RF.

![alt text][LNP]

Finally, we can fit the divS model. This is not trivial. The LNP model is log-concave, meaning that it has a global minimum and no local minima. We only have to run one optimization routine to find it. When fitting the divS model, we initialize with the linear RF and a hacky fit to find a modulator field. We then fit the full modulated poisson model iteratively, conditioning on the RF drive and fitting the divisive Field, then condition on the output of the divisive Field and fit the RF... continue until it converges to something.

This can actually capture the model alright, but I think with some tweaking, we'll do better. We have to build some code to initialize cleverly and learn the hyperparameters (if we need them).

Here's the actual fit fields:
![alt text][divS]

And here's how it does replicating the behavior of the neuron
![alt text][divSsta]

[thedata]: https://github.com/jcbyts/SCdivNorm/tree/master/figures/simulationSTA.png "The Data"

[LNP]: https://github.com/jcbyts/SCdivNorm/tree/master/figures/lnpSTA.png "LNP fit"


[divS]: https://github.com/jcbyts/SCdivNorm/tree/master/figures/divSfits.png "divS fit"

[divSsta]: https://github.com/jcbyts/SCdivNorm/tree/master/figures/divSsta.png "divS STA"
