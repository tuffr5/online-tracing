function new_tracker = reinitialize_tracker(img, bb, tracker)
new_tracker=create_xb_tracker(img, bb);

%REINITIALIZE_HIST_TRACKER Summary of this function goes here
if tracker.use_segmentation
    % retain color hist_gram the main branch hist
    new_tracker.hist_fg = (1-tracker.retain_rate)*new_tracker.hist_fg + tracker.retain_rate*tracker.src_hist_fg;
    new_tracker.hist_bg = (1-tracker.retain_rate)*new_tracker.hist_bg + tracker.retain_rate*tracker.src_hist_bg;
end

% add main branch info
new_tracker.bb6=tracker.bb6;
new_tracker.main_bb6=tracker.main_bb6;
new_tracker.main_c=tracker.main_c;
new_tracker.prev_plane=tracker.prev_plane;
end

