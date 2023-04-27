#!/bin/bash


system="ASSEMBLY.pdb"
topol="topol.top"
index="index.ndx"

gmx grompp -f step6.0_minimization.mdp -o step6.0_minimization.tpr -c ${system} -r ${system} -p ${topol} -n ${index}
gmx mdrun -v -deffnm step6.0_minimization


# Equilibration
cnt=1
cntmax=7

for ((i=${cnt}; i<${cntmax}+1; i++)); do
	if [ $i == "1" ]; then
		gmx grompp -f step6.${i}_equilibration.mdp -o step6.${i}_equilibration.tpr -c step6.$[${i}-1]_minimization.gro -r ${system} -p ${topol} -n ${index}
		gmx mdrun -v -deffnm step6.${i}_equilibration -nt 12
	else
		gmx grompp -f step6.${i}_equilibration.mdp -o step6.${i}_equilibration.tpr -c step6.$[${i}-1]_equilibration.gro -r ${system} -p ${topol} -n ${index}
		gmx mdrun -v -deffnm step6.${i}_equilibration -nt 12
	fi
done
