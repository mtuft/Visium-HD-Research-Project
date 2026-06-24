# modules/queue_module.R

# Protocol: the Shiny process writes a job-manifest RDS to queue_dir/; the worker (worker.R) picks up "queued" jobs, processes them, writes result RDS files and a
# one-line .status file under results_dir/<session_id>/; queueServer() polls every 3s and loads results into the relevant reactiveVals when a job completes.

# Submit a job: writes a manifest RDS (datalist flat-merged with job_id/type/status/timestamp) to the queue and returns the job_id used to poll status.
submitVisiumJob <- function(queue_dir, type, datalist) {
  # [Redacted — manifest construction + write]
}

# Read the worker's one-line status string for a job (NULL if not yet written).
read_worker_status <- function(results_dir, session_id, job_id) {
  # [Redacted]
}

# Queue server: polls on a timer; on "complete" dispatches to the per-job-type completion handler 
# (loads pca/clustered/azimuth/markers/integrated/rctd/subcluster/DE result, updates reactiveVals, advances the active tab, logs);
# on "ERROR:" surfaces the message; renders the job status card; disables run buttons while a job is active; and supports user cancellation.
queueServer <- function(input, output, session,
                        queue_dir, results_dir, session_id,
                        active_job, pca_obj, filtered_obj,
                        markers_data, selection_markers, pipeline_log,
                        update_filtered_obj  = NULL,
                        update_ms_integrated = NULL,
                        ms_de_results        = NULL) {
  # [Redacted — polling, completion dispatch, status card, cancellation]
}

# UI placeholder for the job status card.
queueStatusUI <- function() { uiOutput("queue_status_panel") }
