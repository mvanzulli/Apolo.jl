###################
# Write VTK tests #
###################


start_img = (rand(INTERVAL_START), rand(INTERVAL_START), rand(INTERVAL_START))
length_img = (rand(INTERVAL_LENGTH), rand(INTERVAL_LENGTH), rand(INTERVAL_LENGTH))
offset_img = (0.2, 0.1, 0.3)
finish_img = start_img .+ length_img
