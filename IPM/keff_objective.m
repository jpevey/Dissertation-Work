function [f, g] = keff_objective(variables)
%%% This function takes 121 beta values between 0 - 1, and runs a tsunami
%%% case. When this case is complete the keff (f) and the
%%% derivatives [g] are extracted

pyversion
python_functions = fileparts(which('py_funct.py'))
if count(py.sys.path, python_functions) == 0
    insert(py.sys.path,int32(0), python_functions);
end

py.importlib.import_module('py_funct')

%%% Building the array to pass to python 
python_array = [variables(1),
variables(2),
variables(3),
variables(4),
variables(5),
variables(6),
variables(7),
variables(8),
variables(9),
variables(10),
variables(11),
variables(12),
variables(13),
variables(14),
variables(15),
variables(16),
variables(17),
variables(18),
variables(19),
variables(20),
variables(21),
variables(22),
variables(23),
variables(24),
variables(25),
variables(26),
variables(27),
variables(28),
variables(29),
variables(30),
variables(31),
variables(32),
variables(33),
variables(34),
variables(35),
variables(36),
variables(37),
variables(38),
variables(39),
variables(40),
variables(41),
variables(42),
variables(43),
variables(44),
variables(45),
variables(46),
variables(47),
variables(48),
variables(49),
variables(50),
variables(51),
variables(52),
variables(53),
variables(54),
variables(55),
variables(56),
variables(57),
variables(58),
variables(59),
variables(60),
variables(61),
variables(62),
variables(63),
variables(64),
variables(65),
variables(66),
variables(67),
variables(68),
variables(69),
variables(70),
variables(71),
variables(72),
variables(73),
variables(74),
variables(75),
variables(76),
variables(77),
variables(78),
variables(79),
variables(80),
variables(81),
variables(82),
variables(83),
variables(84),
variables(85),
variables(86),
variables(87),
variables(88),
variables(89),
variables(90),
variables(91),
variables(92),
variables(93),
variables(94),
variables(95),
variables(96),
variables(97),
variables(98),
variables(99),
variables(100),
variables(101),
variables(102),
variables(103),
variables(104),
variables(105),
variables(106),
variables(107),
variables(108),
variables(109),
variables(110),
variables(111),
variables(112),
variables(113),
variables(114),
variables(115),
variables(116),
variables(117),
variables(118),
variables(119),
variables(120),
variables(121)];

%%% This line I call the function evaluate_1d_cyl in the py_funct.py file 
%%% and passes it the list of beta variables.
output = py.py_funct.evaluate_1d_cyl(python_array);

%%% This turns the outuput returned by python into what matlab expects, a
%%% single value and a list of dk/dbeta values.
f = output{1};
g = [ output{2}{1}
output{2}{2}
output{2}{3}
output{2}{4}
output{2}{5}
output{2}{6}
output{2}{7}
output{2}{8}
output{2}{9}
output{2}{10}
output{2}{11}
output{2}{12}
output{2}{13}
output{2}{14}
output{2}{15}
output{2}{16}
output{2}{17}
output{2}{18}
output{2}{19}
output{2}{20}
output{2}{21}
output{2}{22}
output{2}{23}
output{2}{24}
output{2}{25}
output{2}{26}
output{2}{27}
output{2}{28}
output{2}{29}
output{2}{30}
output{2}{31}
output{2}{32}
output{2}{33}
output{2}{34}
output{2}{35}
output{2}{36}
output{2}{37}
output{2}{38}
output{2}{39}
output{2}{40}
output{2}{41}
output{2}{42}
output{2}{43}
output{2}{44}
output{2}{45}
output{2}{46}
output{2}{47}
output{2}{48}
output{2}{49}
output{2}{50}
output{2}{51}
output{2}{52}
output{2}{53}
output{2}{54}
output{2}{55}
output{2}{56}
output{2}{57}
output{2}{58}
output{2}{59}
output{2}{60}
output{2}{61}
output{2}{62}
output{2}{63}
output{2}{64}
output{2}{65}
output{2}{66}
output{2}{67}
output{2}{68}
output{2}{69}
output{2}{70}
output{2}{71}
output{2}{72}
output{2}{73}
output{2}{74}
output{2}{75}
output{2}{76}
output{2}{77}
output{2}{78}
output{2}{79}
output{2}{80}
output{2}{81}
output{2}{82}
output{2}{83}
output{2}{84}
output{2}{85}
output{2}{86}
output{2}{87}
output{2}{88}
output{2}{89}
output{2}{90}
output{2}{91}
output{2}{92}
output{2}{93}
output{2}{94}
output{2}{95}
output{2}{96}
output{2}{97}
output{2}{98}
output{2}{99}
output{2}{100}
output{2}{101}
output{2}{102}
output{2}{103}
output{2}{104}
output{2}{105}
output{2}{106}
output{2}{107}
output{2}{108}
output{2}{109}
output{2}{110}
output{2}{111}
output{2}{112}
output{2}{113}
output{2}{114}
output{2}{115}
output{2}{116}
output{2}{117}
output{2}{118}
output{2}{119}
output{2}{120}
output{2}{121}];
