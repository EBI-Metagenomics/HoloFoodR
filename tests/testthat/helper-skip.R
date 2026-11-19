# Skip tests that require HoloFood database if it is not responding. The
# server is occasionally down, which would otherwise cause false alarms in
# builds.
skip_if_holofood_unavailable <- function(){
    skip_if_offline(host = "holofooddata.org")
    if( !isHoloFoodAvailable() ){
        skip("HoloFood database is not responding.")
    }
}
