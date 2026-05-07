#' Set options for boaR package
#'
#' @param pbStyle An integer from 1 to 3 specifying the progress bar style.
#'   Default is NULL, which will use style 3.
#'   See ?utils::txtProgressBar for details.
#'   When running on an HPC cluster, it is recommended to set pbStyle to 1 for better output.
#'
#' @export
#' @examples
#' # Set the progress bar style to 1 for the current session
#' # set_boaR_options(pbStyle = 1)

set_boaR_options <- function(pbStyle = NULL) {
	if (!is.null(pbStyle)) {
		if (!(pbStyle %in% 1:3)) {
			stop("pbStyle must be an integer between 1 and 3.", call. = FALSE)
		}
		options(pbStyle = pbStyle)
	}

	invisible()
}
