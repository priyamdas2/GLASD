function y = eval_count_store(val)
    global eval_count eval_history
    if isvector(val)
        eval_count = eval_count + 1;
        eval_history(end+1) = val;
    else
        eval_count = eval_count + size(val,1);
        eval_history = [eval_history; val(:)];
    end
    y = val;
end