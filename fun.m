function fitness = fun(Parameters, Ft, target)

predict = energy_model(Parameters, Ft);

fitness = abs(predict-target);

end