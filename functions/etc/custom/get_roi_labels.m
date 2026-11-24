function roi_labels = get_roi_labels(p)

% Extract labels
temp = load(p.sourcemodel.grid.ts.ve(p.subject(1), p.session(1), p.run(1),p.task(1)));
roi_labels = temp.pos_name;

end