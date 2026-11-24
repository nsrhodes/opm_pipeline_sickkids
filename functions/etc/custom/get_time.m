function time = get_time(p)

% Extract time
temp = load(p.sourcemodel.grid.ts.ve(p.subject(1), p.session(1), p.run(1),p.task(1)));
time = temp.time';

end