
disp('1. Two Loop Network')
disp('2. Hanoi Network')
disp('3. Newyork City Network')
disp('4. GoYang Network')
disp('5. BakRyun Network')
n = input('Enter a number: ');

switch n
    case 1
        d=epanet('twoloop.inp');
        d.plot;
    case 2
        d=epanet('hanoi.inp');
        d.plot;
    case 3
        d=epanet('newyork.inp');
        d.plot;
    case 4
        d=epanet('goyang.inp');
        d.plot;
    case 5 
        d=epanet('bakryun.inp');
        d.plot;
    otherwise
        disp('Enter Correct Number')
end