// Vapour saturation pressure and change
{
    dimensionedScalar a0("temp", dimPressure, 610.5);
    dimensionedScalar a1("temp", dimless, 17.269);
    dimensionedScalar a2("temp", dimTemperature, 4717.03);
    dimensionedScalar a3("temp", dimTemperature, 35.85);
    ps = a0*pow(2.718,((a1*T - a2) / (T - a3)));
}

{
    dimensionedScalar a0("temp", dimPressure/dimTemperature, 610.5);
    dimensionedScalar a1("temp", dimTemperature, 4717.03);
    dimensionedScalar a2("temp", dimTemperature, 35.85);
    dimensionedScalar a3("temp", dimTemperature, 17.2672094763233);
    dimensionedScalar a4("temp", dimTemperature, 0.999896315728952);
    dpsDT = a0*pow(2.718, (17.269*T - a1)/(T - a2))
               *(a3/(T - a2) 
                - a4*(17.269*T - a1)/pow(T - a2,2)
                 );
}
