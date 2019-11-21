# Install headers
file(GLOB relnav_headers "*.h")
install(FILES ${relnav_headers} DESTINATION include/gtsam/relnav)

# Build tests
add_subdirectory(tests)
