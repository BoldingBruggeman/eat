# Copyright (C) 2021 Bolding & Bruggeman

#KBadd_subdirectory(seamless)
#KBadd_custom_target(test_eat_seamless DEPENDS eat_seamless_filter DEPENDS eat_seamless_model DEPENDS eat_seamless_obs)

add_subdirectory(2d)
add_custom_target(test_eat_2d DEPENDS eat_pdaf_filter DEPENDS eat_2d_obs DEPENDS eat_2d_model)

#KBadd_custom_target(test_eat_all DEPENDS test_eat_2d DEPENDS test_eat_seamless)
add_custom_target(test_eat_all DEPENDS test_eat_2d)
