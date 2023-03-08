% output from RoSE
%     [  7  ]
%     [ 2 3 ]
%     [6   8]
%     [ 1 4 ]
%     [  5  ]
% input for channel registration
%     [  1  ]
%     [ 8 2 ]
%     [7   3]
%     [ 6 4 ]
%     [  5  ]

function x = RoSE2ChReg(x0)
x{1} = x0{7};
x{2} = x0{3};
x{3} = x0{8};
x{4} = x0{4};
x{5} = x0{5};
x{6} = x0{1};
x{7} = x0{6};
x{8} = x0{2};
end

