import amber_adaptbx

def get_amber_struct_object(params):
  amber_params = params.amber
  sander=False
  if params.amber.md_engine == "sander":
    import sander, sanderles
    ridingH=True
    if getattr(params, "hydrogens", False): # ensemble refinement does not have this
        if params.hydrogens.refine in ['riding', 'Auto']:
          ridingH = True
        elif params.hydrogens.refine in ['individual']:
          ridingH = False
        else:
          raise Sorry("Hydrogens.refine parameter '%s' unknown!"
                      %params.hydrogens.refine)
    amber_structs = amber_adaptbx.sander_structs(
      parm_file_name=amber_params.topology_file_name,
      rst_file_name=amber_params.coordinate_file_name,
      ridingH=ridingH,
      )
    if amber_structs.is_LES:
      sanderles.setup(amber_structs.parm,
             amber_structs.rst.coordinates,
             amber_structs.rst.box,
             amber_structs.inp)
    else:
      sander.setup(amber_structs.parm,
             amber_structs.rst.coordinates,
             amber_structs.rst.box,
             amber_structs.inp)
    sander=True
  elif params.amber.md_engine == "mdgx":
    amber_structs = amber_adaptbx.mdgx_structs(
      parm_file_name=amber_params.topology_file_name,
      rst_file_name=amber_params.coordinate_file_name)
  else:
    raise Sorry("Unsupported md_engine %s" \
          %params.amber.md_engine)

  return amber_structs, sander
