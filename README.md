The shared codes and data include necessary materials for reproducing the results of our paper. They are organized according to the order of figs in the paper.

Single-lane ring:

Fig1: 
fig1_a.ipynb: simulation in SUMO without autonomous vehicle (AV)
fig1_b.ipynb: simulation in SUMO with AV

Fig2:
fig2.m: speed profiles and trajectories of all human-driven vehicles and the AV when controlled by the benign model

Fig3:
fig3.m: comparison between probability densities and observed frequencies of the random variables (v_max, d_crit, d_min, T_mix) for trigger selection

Fig4:
fig4.m: speed profiles and trajectories of all human-driven vehicles and the AV when controlled by the malicious model

Fig5:
fig5.m: speed profiles and trajectories of all human-driven vehicles and the AV in congestion attack

Fig6:
fig6.m: speed profiles and trajectories of all human-driven vehicles and the AV for the benign model when the trigger sample is encountered

Fig7:
fig7.ipynb: simulation in SUMO of insurance attack

Trigger_exploration:
trigger_exploration.m: trigger exploration for the insurance attack and the insurance attack


Neural networks: 

Containing the benign and backdoor models for single-lane ring (single_lane_notrigger.sav, single_lane_congestion.sav, single_lane_insurance.sav)



In our paper, the training of benign models are based on the DDPG algorithm and thanks to the ddpg agent building by Chris Yoon on GitHub. 




