#! /bin/bash

cat >params.eff <<EOF
refinement {
  main {
    number_of_macro_cycles=1
  }
  refine {
    strategy = *individual_sites individual_sites_real_space rigid_body \
               *individual_adp group_adp tls *occupancies group_anomalous
    adp {
      individual {
        anisotropic = "not (element H)"
        isotropic = "element H"
      }
    }
  }
  bulk_solvent_and_scale {
    bulk_solvent = True
    anisotropic_scaling = True
    k_sol_b_sol_grid_search = False
    minimization_k_sol_b_sol = False
  }
}
EOF

phenix.refine \
	params.eff \
	4lzt.pdb \
	4lzt-sf-truncated.mtz \
	refinement.input.xray_data.r_free_flags.file_name=md_avg_95_rfree.mtz \
	use_amber=True \
	topology_file_name=4lzt.prmtop \
	coordinate_file_name=4lzt.rst7 \
	md_engine=sander \
	wxc_scale=0.025 \
	prefix=vs_obs_amber \
	--overwrite

