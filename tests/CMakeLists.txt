# Copyright (C) 2021 Bolding & Bruggeman

add_custom_target(all_handlers)
#KBadd_dependencies(all_handlers eat_server eat_model eat_observations)
add_dependencies(all_handlers eat_model eat_observations)

# mpiexec -np 1 eat_server : -np 1 eat_observation :  -np 6 eat_worker

add_test(NAME help COMMAND eat_server -h)
