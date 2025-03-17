
format:
	@Rscript --quiet -e "library(styler); style_dir('R'); style_file('ui.R'); style_file('server.R');"


systemd:
	@footprint config template etc/lfqanalyst.service r_home=${CONDA_PREFIX} port=4447 \
		data-dir=./data-dir

nginx:
	@footprint config template etc/lfqanalyst.conf \
	 	server-name=lfq.plantenergy.org port=4447
# documentation
oxygen:
	@Rscript --quiet -e "roxygen2::roxygenize()"

.PHONY: format oxygen systemd nginx
