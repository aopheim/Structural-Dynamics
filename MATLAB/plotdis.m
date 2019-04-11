function  outputstring=plotdis(Nodes,Elements,locnod,disp,ampl),

% plots the mesh with displacement multiply by ampl

Nel = size(Elements,1);

%plot3(Nodes(:,1),Nodes(:,2),Nodes(:,3),'g*')
for el=1:Nel,
    type   = Elements(el,1);
    hold on
    dof1 = locnod(Elements(el,3),[1 2 3]);
    dof2 = locnod(Elements(el,4),[1 2 3]);
    %add the displacements to the coordinates of the nodes of the element
    node1  = Nodes(Elements(el,3),:) + disp(dof1)'*ampl;
    node2  = Nodes(Elements(el,4),:) + disp(dof2)'*ampl;
    if type==1,
        plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],'k', 'Linewidth', 2);
    elseif type==2,
        plot3([node1(1);node2(1)],[node1(2);node2(2)],[node1(3);node2(3)],'k', 'Linewidth', 2);
    end
    clear dof1 dof2
end
axis equal
axis off
view(-40,15)
%view(0,0)